import xarray as xr
import numpy as np
from numba import float64, guvectorize
import os
import pandas as pd
from pywapor.general.processing_functions import save_ds
from pywapor.enhancers.smooth.plotters import plot_overview, plot_stats, create_video
import tqdm
import glob
from joblib import Parallel, delayed

def open_ts(fps):

    dss = list()
    for fp in fps:
        name = os.path.split(fp)[-1].replace(".nc", "")
        ds = xr.open_dataset(fp, decode_coords = "all")
        ds["sensor"] = xr.ones_like(ds["time"], dtype = int).where(False, name)
        dss.append(ds)
    ds = xr.concat(dss, dim = "time", combine_attrs = "drop_conflicts").sortby("time").transpose("y", "x", "time")
    ds["sensor"] = ds["sensor"].astype("<U7")
    ds.attrs = {}

    attribute = {str(i): sensor_name for i, sensor_name in enumerate(np.unique(ds.sensor))}
    values = np.array(list(attribute.keys()), dtype = int)
    coords = np.array(list(attribute.values()), dtype = str)
    transformer = xr.DataArray(values, dims=["sensor"], coords = {"sensor": coords})
    ds["sensor"] = transformer.sel(sensor = ds.sensor).drop("sensor").assign_attrs(attribute)

    return ds

def filter_ts(x, y, tol = 0.002):
    diff_forward = x.diff("time").assign_coords({"time": x.time[:-1]})
    diff_backward = x.diff("time").assign_coords({"time": x.time[1:]})
    d2x = (diff_forward + diff_backward).reindex_like(x.time, fill_value = tol * 2)
    print(f"--> Removing {(d2x < tol).sum().values} time slices.")
    y = y.where(d2x >= tol, drop = True)
    x = x.where(d2x >= tol, drop = True)
    return x, y

def choose_func(y, lmbdas, fname):

    funcs = {"wt": [wt1, wt2], "cve": [cve1, cve2]}

    y_dims = getattr(y, "ndim", 0)
    lmbdas_dims = getattr(lmbdas, "ndim", 0)
    if y_dims in [2,3] and lmbdas_dims in [1]:
        wt = funcs[fname][1]
    elif y_dims in [2] and lmbdas_dims in [2]:
        raise ValueError
    else:
        wt = funcs[fname][0]

    if y_dims == 3 and lmbdas_dims == 2 and not isinstance(y, xr.DataArray):
        assert y.shape[:2] == lmbdas.shape
    elif y_dims == 3 and lmbdas_dims == 2 and isinstance(y, xr.DataArray):
        assert np.all([True for k, v in lmbdas.sizes.items() if y.sizes[k] == v])

    print(f"Using {wt.__name__}")

    return wt

def second_order_diff_matrix(x):
    X = np.append(x, [np.nan, np.nan])
    # Create x-aware delta matrix. When sample-points are at equal distance,
    # this reduces to [[1, -2, 1, ...], ..., [..., 1, -2, 1]].
    diag_n0 = 2  / ((X[1:-1] - X[:-2])  * (X[2:]   - X[:-2]))
    diag_n1 = -2 / ((X[2:-1] - X[1:-2]) * (X[1:-2] - X[0:-3]))
    diag_n2 = 2  / ((X[2:-2] - X[1:-3]) * (X[2:-2] - X[0:-4]))
    D = np.diag(diag_n0) + np.diag(diag_n1, k = 1) + np.diag(diag_n2, k = 2)
    D = D[:-2,:]
    return D

def whittaker(y, x, lmbdas = 10.):

    x, y, lmbdas, dim_name = assert_dtypes(x, y, lmbdas)

    # Normalize x-coordinates
    x = (x - x.min()) / (x.max() - x.min()) * x.size

    # Remove points that are too close together
    if isinstance(x, xr.DataArray):
        original_time = x[dim_name]
        x, y = filter_ts(x, y, tol = 0.002)

    # Create x-aware delta matrix.
    D = second_order_diff_matrix(x)
    A = np.dot(D.T, D)

    # Choose which Whittaker function to use depending on shapes of y and lmbdas.
    wt = choose_func(y, lmbdas, "wt")

    if isinstance(x, xr.DataArray):

        if wt.__name__ == "wt1":
            icd = [[dim_name], [], []]
            ocd = [[dim_name]]
        elif wt.__name__ == "wt2":
            
            icd = [[dim_name], [], ["lmbda"]]
            ocd = [["lmbda", dim_name]]

        # Make sure lmbdas is chunked similar to y.
        if not isinstance(y.chunk, type(None)):
            lmbdas = lmbdas.chunk({k: v for k,v in y.chunksizes.items() if k in lmbdas.dims})

        # Apply whittaker smoothing along axis.
        z = xr.apply_ufunc(
            wt, y, A, lmbdas,
            input_core_dims=icd,
            output_core_dims=ocd,
            dask = "allowed",
            )

        z = z.reindex_like(original_time, fill_value = np.nan)
    else:
        z = wt(y, A, lmbdas)

    return z

@guvectorize([(float64[:], float64[:,::1], float64[:], float64[:])], '(n),(n,n),()->(n)')
def wt1(Y, A, lmbda, Z):
    # Create weights.
    w = np.isfinite(Y, np.zeros_like(Y))
    W = np.diag(w)
    # Set np.nan in y to 0.
    Y = np.where(w == 0, 0, Y)
    Z[:] = np.linalg.solve(W + lmbda * A, np.dot(W,Y))

@guvectorize([(float64[:], float64[:,::1], float64[:], float64[:,:])], '(n),(n,n),(k)->(k,n)')
def wt2(Y, A, lmbda, Z):
    # Create weights.
    w = np.isfinite(Y, np.zeros_like(Y))
    W = np.diag(w)
    # Set np.nan in y to 0.
    Y = np.where(w == 0, 0, Y)
    for i, lmbda in enumerate(lmbda):
        Z[i, :] = np.linalg.solve(W + lmbda * A, np.dot(W,Y))

def cross_val_lmbda(y, x, lmbdas = np.logspace(0, 3, 4)):

    x, y, lmbdas, dim_name = assert_dtypes(x, y, lmbdas)

    # Normalize x-coordinates
    x = (x - x.min()) / (x.max() - x.min()) * x.size

    # Remove points that are too close together
    if isinstance(x, xr.DataArray):
        x, y = filter_ts(x, y, tol = 0.002)

    # Create x-aware delta matrix.
    D = second_order_diff_matrix(x)
    A = np.dot(D.T, D)

    # Choose which function to use depending on shapes of y and lmbdas.
    cve = choose_func(y, lmbdas, "cve")

    if isinstance(x, xr.DataArray):

        # Determine dimension names.
        if cve.__name__ == "cve1":
            icd = [[], [dim_name], []]
            ocd = [[]]
        elif cve.__name__ == "cve2":
            icd = [["lmbda"], [dim_name], []]
            ocd = [["lmbda"]]

        # Make sure lmbdas is chunked similar to y.
        if not isinstance(y.chunk, type(None)):
            lmbdas = lmbdas.chunk({k: v for k,v in y.chunksizes.items() if k in lmbdas.dims})

        # Calculate the cross validation standard error for each lambda using
        # the hat matrix.
        cves = xr.apply_ufunc(
            cve, lmbdas, y, A,
            input_core_dims=icd,
            output_core_dims=ocd,
            dask = "allowed")

        # Select the lambda for which the error is smallest.
        if "lmbda" in cves.dims:
            lmbda_sel = cves.idxmin(dim = "lmbda")
        else:
            cves.assign_coords({"lmbda": lmbdas})
            lmbda_sel = lmbdas
    else:
        cves = cve(lmbdas, y, A)
        if np.isscalar(lmbdas):
            lmbda_sel = lmbdas
        elif lmbdas.ndim == 1:
            idx = np.argmin(cves, axis = -1)
            lmbda_sel = lmbdas[idx]
        else:
            lmbda_sel = lmbdas
        
    return lmbda_sel, cves

@guvectorize([(float64, float64[:], float64[:,::1], float64[:])], '(),(m),(m,m)->()')
def cve1(lmbda, Y, A, cves):
    w = np.isfinite(Y, np.zeros_like(Y))
    W = np.diag(w)
    Y_ = np.where(w == 0, 0, Y)
    H = np.linalg.solve(W + lmbda * A, W) # Eq. 13
    y_hat = np.dot(H, Y_) # Eq. 10
    hii = np.diag(H)
    cves[:] = np.sqrt(np.nanmean(((Y - y_hat)/(1- hii))**2)) # Eq. 9 + 11

@guvectorize([(float64[:], float64[:], float64[:,::1], float64[:])], '(k),(m),(m,m)->(k)')
def cve2(lmbdas, Y, A, cves):
    w = np.isfinite(Y, np.zeros_like(Y))
    W = np.diag(w)
    Y_ = np.where(w == 0, 0, Y)
    for i, lmbda in enumerate(lmbdas):
        H = np.linalg.solve(W + lmbda * A, W) # Eq. 13
        y_hat = np.dot(H, Y_) # Eq. 10
        hii = np.diag(H)
        cves[i,...] = np.sqrt(np.nanmean(((Y - y_hat)/(1- hii))**2)) # Eq. 9 + 11

def assert_dtypes(x, y, lmbdas):

    # Check x and y.
    assert x.ndim == 1
    if isinstance(x, xr.DataArray) and isinstance(y, xr.DataArray):
        dim_name = x.dims[0]
        assert dim_name in y.dims
    elif isinstance(x, np.ndarray) and isinstance(y, xr.DataArray):
        dim_names = [k for k,v in y.sizes.items() if v == x.size]
        if len(dim_names) != 1:
            raise ValueError
        else:
            dim_name = dim_names[0]
            x = xr.DataArray(x, dims = [dim_name], coords = y.dim_name)
    elif isinstance(x, xr.DataArray) and isinstance(y, np.ndarray):
        x = x.values
        dim_name = None
        assert x.size == y.shape[-1]
    elif isinstance(x, np.ndarray) and isinstance(y, np.ndarray):
        assert x.size == y.shape[-1]
        dim_name = None
    else:
        raise TypeError

    # Check lmbdas.
    assert lmbdas.ndim <= 2
    if isinstance(x, xr.DataArray) and (isinstance(lmbdas, np.ndarray) or np.isscalar(lmbdas)):
        if not np.isscalar(lmbdas):
            assert lmbdas.ndim <= 1
            if lmbdas.ndim == 0:
                lmbdas = float(lmbdas)
            else:
                lmbdas = xr.DataArray(lmbdas, dims = ["lmbda"], coords = {"lmbda": lmbdas})
        # else:
        lmbdas = xr.DataArray(lmbdas)
    elif isinstance(x, xr.DataArray) and isinstance(lmbdas, xr.DataArray):
        ...
    elif isinstance(x, np.ndarray) and (isinstance(lmbdas, np.ndarray) or np.isscalar(lmbdas)):
        if lmbdas.ndim == 0:
            lmbdas = float(lmbdas)
        elif lmbdas.ndim == 2 and y.ndim == 3:
            assert y.shape[:-1] == lmbdas.shape
    elif isinstance(x, np.ndarray) and isinstance(lmbdas, xr.DataArray):
        lmbdas = lmbdas.values
        if lmbdas.ndim == 0:
            lmbdas = float(lmbdas)
    else:
        raise TypeError
    
    return x, y, lmbdas, dim_name


def whittaker_smoothing(ds, var, lmbdas = 100., out_fh = None, xdim = "time", new_x = None, export_all = False):

    if isinstance(out_fh, type(None)):
        folder = os.path.split(ds.encoding["source"])[0]
        out_fh = os.path.join(folder, f"smoothed_{var}.nc")
    else:
        folder = os.path.split(out_fh)[0]

    if not os.path.exists(folder):
        os.makedirs(folder)

    while os.path.isfile(out_fh):
        name, ext = os.path.splitext(out_fh)
        out_fh = name + "_" + ext

    # Add new x values.
    if not isinstance(new_x, type(None)) and getattr(xdim, '__len__', lambda: 0)() > 0:
        ds = xr.merge([ds, xr.Dataset({xdim: new_x})]).sortby(xdim).chunk({xdim: -1})
        if "sensor" in ds.data_vars:
            sensor_id = np.max(np.unique(ds["sensor"])) + 1
            ds["sensor"] = ds["sensor"].fillna(sensor_id).assign_attrs({str(sensor_id): "Interp."})

    ds = ds.assign_coords({"lmbda": lmbdas})

    # Only do this when more then one lmbda is provided.
    if getattr(lmbdas, 'size', 1) > 1 or export_all:
        ds["lmbda_sel"], ds["cves"] = cross_val_lmbda(ds[var], ds[xdim], lmbdas = lmbdas)
        lmbdas = ds["lmbda_sel"]

    # Smooth the data.
    ds[f"{var}_smoothed"] = whittaker(ds[var], ds[xdim], lmbdas = lmbdas)

    if not export_all:
        ds = ds[[f"{var}_smoothed"]].rename({f"{var}_smoothed": var})
        if not isinstance(new_x, type(None)) and getattr(xdim, '__len__', lambda: 0)() > 0:
            ds = ds.sel({xdim:  new_x})

    ds = save_ds(ds, out_fh, encoding = "initiate", scheduler = "synchronous")

    return ds

def make_plots(ds, folder, points, xdim, new_x, parallel = True):

    if not os.path.exists(folder):
        os.makedirs(folder)

    if not parallel:
        for i in tqdm.tqdm(np.arange(0, ds[xdim].size, 1)[~np.isin(ds[xdim], new_x)]):
            plot_overview(ds, points, i, folder)
    else:
        _ = Parallel(n_jobs=4)(delayed(plot_overview)(ds, points, i, folder) for i in np.arange(0, ds[xdim].size, 1)[~np.isin(ds[xdim], new_x)])
    files = np.sort(glob.glob(os.path.join(folder, "[0-9]" * 6 + ".png")))
    video_fh = os.path.join(folder, "timeseries.mp4")
    create_video(files, video_fh)
    plot_stats(ds, folder)

def create_points(ds, n, offset = 1):
    lons = ds.x.isel(x=np.linspace(0+offset, ds.x.size-(1+offset), n, dtype=int)).values
    lats = ds.y.isel(y=np.linspace(0+offset, ds.y.size-(1+offset), n, dtype=int)).values
    ys, xs = np.meshgrid(lats, lons)
    ys = ys.flatten().tolist()
    xs = xs.flatten().tolist()
    names = [f"P{i:>02}" for i in range(1, len(xs)+1)]
    return (xs, ys, names)

if __name__ == "__main__":

    ts_fp =  r"/Users/hmcoerver/Local/input_test_series.nc"

    ls7_fp = r"/Users/hmcoerver/Local/long_timeseries/fayoum/LANDSAT/LE07_SR.nc"
    ls8_fp = r"/Users/hmcoerver/Local/long_timeseries/fayoum/LANDSAT/LC08_SR.nc"
    ls9_fp = r"/Users/hmcoerver/Local/long_timeseries/fayoum/LANDSAT/LC09_SR.nc"
    mod13_fp = r"/Users/hmcoerver/Local/long_timeseries/fayoum/MODIS/MOD13Q1.061.nc"
    myd13_fp = r"/Users/hmcoerver/Local/long_timeseries/fayoum/MODIS/MYD13Q1.061.nc"

    fps = [
        ls7_fp,
        ls8_fp,
        ls9_fp,
        # mod13_fp,
        # myd13_fp,
    ]
    
#     # ds = xr.open_dataset(ts_fp, decode_coords="all")
#     # ds = open_ts(fps).sel(x = slice(30.812, 30.825), y = slice(29.450, 29.44))#.isel(x = slice(200, 225), y = slice(200, 220))
#     ds = open_ts(fps).isel(x = slice(300, 320), y = slice(500, 520))#.isel(x = slice(200, 225), y = slice(200, 220))
#     # ds = open_ts(fps).isel(x = 310, y = 505)#.isel(x = slice(200, 225), y = slice(200, 220))
    ds = open_ts(fps)
    ds = ds.chunk({"time": -1, "x": -1, "y": -1})

    folder = r"/Users/hmcoerver/Local/test_02" #
    out_fh = os.path.join(folder, "test_03.nc") #
#     # lats = [29.46376753]
#     # lons = [30.78233174]
#     # names = ["P1"]
#     # points = (lons, lats, names)
#     points = create_points(ds, 2, offset = 1)
#     # points = None
#     lmbdas = np.logspace(1,4,4)
#     new_x = pd.date_range(np.datetime_as_string(ds.time[0], "D") + " 12:00", 
#                           np.datetime_as_string(ds.time[-1], "D") + " 12:00", freq="D")
#     # new_x = None
#     var = "ndvi"
#     xdim = "time"
#     export_all = True

    # out = whittaker_smoothing(ds, var, lmbdas = lmbdas, out_fh = out_fh, new_x = new_x, points = points, export_all=export_all)

    # make_plots(out, os.path.join(folder, "graphs"), points, xdim, new_x, parallel = True)