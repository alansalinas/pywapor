import xarray as xr
import numpy as np
from numba import float64, guvectorize
import os
import pandas as pd
from pywapor.general.processing_functions import save_ds
from pywapor.enhancers.smooth.plotters import plot_overview

def open_ts(fps):
    dss = list()
    for fp in fps:
        name = os.path.split(fp)[-1].replace(".nc", "")
        ds = xr.open_dataset(fp, decode_coords = "all")
        ds["sensor"] = xr.ones_like(ds["time"], dtype = int).where(False, name)
        dss.append(ds)
    ds = xr.concat(dss, dim = "time", combine_attrs = "drop_conflicts").sortby("time").transpose("y", "x", "time")
    ds.attrs = {}
    return ds

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

def second_order_diff_matrix(x, normalize = True):
    if normalize:
        X = (x - x.min()) / (x.max() - x.min()) * x.size
    else:
        X = x[:]
    X = np.append(X, [np.nan, np.nan])
    # Create x-aware delta matrix. When sample-points are at equal distance,
    # this reduces to [[1, -2, 1, ...], ..., [..., 1, -2, 1]].
    diag_n0 = 2  / ((X[1:-1] - X[:-2])  * (X[2:]   - X[:-2]))
    diag_n1 = -2 / ((X[2:-1] - X[1:-2]) * (X[1:-2] - X[0:-3]))
    diag_n2 = 2  / ((X[2:-2] - X[1:-3]) * (X[2:-2] - X[0:-4]))
    D = np.diag(diag_n0) + np.diag(diag_n1, k = 1) + np.diag(diag_n2, k = 2)
    D = D[:-2,:]
    return D

def whittaker(y, x, lmbdas = 10., axis = -1):

    assert y.shape[axis] == x.size

    if isinstance(y, xr.DataArray) and not isinstance(lmbdas, xr.DataArray):
        lmbdas = xr.DataArray(lmbdas)

    if x.dtype != int:
        x = x.astype(int)

    # Make sure correct axis is smoothed.
    if isinstance(axis, int) and axis != -1:
        y = np.moveaxis(y, axis, -1)
    elif isinstance(axis, str):
        assert isinstance(y, xr.DataArray)
    
    # Determine y core dimension name.
    if isinstance(axis, int):
        dummy_names = [f"dim_{i}" for i in range(getattr(y, "ndim", 0))]
        dim_name = getattr(y, "dims", dummy_names)[axis]
    elif isinstance(axis, str):
        dim_name = axis

    # Choose which Whittaker function to use depending on shapes of y and lmbdas.
    wt = choose_func(y, lmbdas, "wt")

    # Create x-aware delta matrix.
    D = second_order_diff_matrix(x)
    A = np.dot(D.T, D)

    if wt.__name__ == "wt1":
        icd = [[dim_name], [], []]
        ocd = [[dim_name]]
    elif wt.__name__ == "wt2":
        lmbda_dim_name = getattr(lmbdas, "dims", ["lmbda"])[0]
        icd = [[dim_name], [], [lmbda_dim_name]]
        ocd = [[lmbda_dim_name, dim_name]]

    # Apply whittaker smoothing along axis.
    z = xr.apply_ufunc(
        wt, y, A, lmbdas,
        input_core_dims=icd,
        output_core_dims=ocd,
        # vectorize = True,
        dask = "allowed",
        # output_dtypes=[y.dtype],
        )

    # Put axes back in original order.
    if isinstance(axis, int) and axis != -1:
        z = np.moveaxis(z, -1, axis)

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

    if isinstance(y, xr.DataArray) and not isinstance(lmbdas, xr.DataArray):
        lmbdas = xr.DataArray(lmbdas, dims = ["lmbda"], coords = {"lmbda": lmbdas})

    if x.dtype != int:
        x = x.astype(int)

    # Create x-aware delta matrix.
    D = second_order_diff_matrix(x)
    A = np.dot(D.T, D)

    # Choose which function to use depending on shapes of y and lmbdas.
    cve = choose_func(y, lmbdas, "cve")

    # Determine dimension names.
    dummy_names = [f"dim_{i}" for i in range(getattr(lmbdas, "ndim", 0))]
    if cve.__name__ == "cve1":
        lmbda_dim_names = getattr(lmbdas, "dims", dummy_names)
        lmbda_dim_name = lmbda_dim_names[0] if 0 < len(lmbda_dim_names) < 2 else None
        icd = [[], ["time"], []]
        ocd = [[]]
    elif cve.__name__ == "cve2":
        lmbda_dim_name = getattr(lmbdas, "dims", dummy_names)[0]
        icd = [[lmbda_dim_name], ["time"], []]
        ocd = [[lmbda_dim_name]]

    # Calculate the cross validation standard error for each lambda using
    # the hat matrix.
    cves = xr.apply_ufunc(
        cve, lmbdas, y, A,
        input_core_dims=icd,
        output_core_dims=ocd,
        # vectorize = True,
        dask = "allowed",
        # output_dtypes=[ds.ndvi.dtype],
        )

    # Select the lambda for which the error is smallest.
    if not isinstance(lmbda_dim_name, type(None)):
        if isinstance(cves, xr.DataArray):
            lmbda_sel = cves.idxmin(dim = lmbda_dim_name)
        else:
            lmbda_sel = lmbdas[np.argmin(cves, axis = -1)]
    else: 
        lmbda_sel = None
        if isinstance(cves, xr.DataArray):
            cves = cves.assign_coords({"lmbda": lmbdas})

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

def whittaker_smoothing(ds, var, lmbdas = 100., out_fh = None, xdim = "time", new_x = None, points = None):

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
        ds = xr.merge([ds, xr.Dataset({xdim: new_x})]).sortby(xdim)
        if "sensor" in ds.data_vars:
            ds["sensor"] = ds["sensor"].fillna("Interp.")

    # Only do this when more then one lmbda is provided.
    if getattr(lmbdas, '__len__', lambda: 1)() > 1:
        ds["lmbda_sel"], ds["cves"] = cross_val_lmbda(ds[var], ds[xdim], lmbdas = lmbdas)
        lmbdas = ds["lmbda_sel"]

    # Smooth the data.
    ds[f"{var}_smoothed"] = whittaker(ds[var], ds[xdim], lmbdas = lmbdas)

    if not isinstance(points, type(None)):
        ds = ds.compute(scheduler = 'synchronous')
        for i in np.arange(0, ds[xdim].size, 1)[~np.isin(ds[xdim], new_x)]:
            plot_overview(ds, points, i, folder)

    if not isinstance(new_x, type(None)) and getattr(xdim, '__len__', lambda: 0)() > 0:
        ds = ds.sel({xdim:  new_x})
        
    ds = ds[[f"{var}_smoothed"]].rename({f"{var}_smoothed": var})
    
    ds = save_ds(ds, out_fh, encoding = "initiate", scheduler = "synchronous")

    return ds

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
    
    # ds = xr.open_dataset(ts_fp, decode_coords="all")
    # ds = open_ts(fps).sel(x = slice(30.812, 30.825), y = slice(29.450, 29.44))#.isel(x = slice(200, 225), y = slice(200, 220))
    ds = open_ts(fps).isel(x = slice(300, 310), y = slice(500, 510))#.isel(x = slice(200, 225), y = slice(200, 220))
    # ds = open_ts(fps).isel(x = 310, y = 505)#.isel(x = slice(200, 225), y = slice(200, 220))
    # ds = open_ts(fps)
    ds = ds.chunk({"time": -1, "x":-1, "y":-1})

    folder = r"/Users/hmcoerver/Local/test_01"#
    out_fh = os.path.join(folder, "test_03.nc")#
    make_plots = False
    lats = [29.462, 29.444]
    lons = [30.784, 30.787]
    names = ["P1", "P2"]
    points = (lons, lats, names)
    lmbdas = 100.
    new_x = pd.date_range(np.datetime_as_string(ds.time[0], "D") + " 12:00", 
                          np.datetime_as_string(ds.time[-1], "D") + " 12:00", freq="D")
    var = "ndvi"
    xdim = "time"

    # out = whittaker_smoothing(ds, var, new_x = new_x, points = points)


