import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay, cKDTree
from dask.diagnostics import ProgressBar

def create_grid(input_ds, dx, dy):

    # The left and lower bounds will snap to a `1/precision`-th degree 
    # below smallest coordinate in input_ds.
    precision = 0.1

    # Open input dataset if necessary.
    if isinstance(input_ds, str):
        input_ds = xr.open_dataset(input_ds, chunks = "auto")

    # Determine box that is just larger then domain in input_ds.
    xmin = np.floor(input_ds.x.min().values * precision**-1) / precision**-1
    xmax = np.ceil(input_ds.x.max().values * precision**-1) / precision**-1
    ymin = np.floor(input_ds.y.min().values * precision**-1) / precision**-1
    ymax = np.ceil(input_ds.y.max().values * precision**-1) / precision**-1

    # Determine coordinates.
    gridx = np.arange(xmin, xmax + dx, dx)
    gridy = np.arange(ymin, ymax + dy, dy)

    # Wrap coordinates in xr.Dataset.
    grid_ds = xr.Dataset(None, coords = {"y": gridy, "x": gridx})
    
    # Add relevant variables.
    for var in input_ds.data_vars:
        dims = [(x, input_ds[x]) for x in input_ds[var].dims if x not in ("ny", "nx")]
        dims += [("y", grid_ds.y), ("x", grid_ds.x)]
        grid_ds[var] = xr.DataArray(coords = {k:v for k,v, in dims}, dims = [x[0] for x in dims])

    return grid_ds

def regrid(grid_ds, input_ds, max_px_dist = 10):

    # Create output dataset
    output_ds = grid_ds.stack({"grid_pixel": ("y", "x")})
    no_pixel = output_ds.grid_pixel.size

    # Create intermediate dataset without multi-index
    grid_adj_ds = output_ds.assign_coords({"i": ("grid_pixel", range(no_pixel))}).set_index(grid_pixel="i")
    grid_adj_ds = grid_adj_ds.drop(["x", "y"]).assign({"x": ("grid_pixel", grid_adj_ds.x.values),
                                                       "y": ("grid_pixel", grid_adj_ds.y.values)})

    # Determine bounding box of current chunk.
    bb = [grid_adj_ds.x.min(), grid_adj_ds.y.min(), grid_adj_ds.x.max(), grid_adj_ds.y.max()]

    # Determine pixel size.
    dx = np.median(np.unique(np.diff(grid_ds.x)))
    dy = np.median(np.unique(np.diff(grid_ds.y)))

    # Open input dataset if necessary.
    if isinstance(input_ds, str):
        input_ds = xr.open_dataset(input_ds, chunks = "auto")

    # Filter out irrelevant input data for current chunk.
    ymask = input_ds.y.where((input_ds.y >= bb[1] - dy) &
                             (input_ds.y <= bb[3] + dy), drop = False)
    xmask = input_ds.x.where((input_ds.x >= bb[0] - dx) & 
                             (input_ds.x <= bb[2] + dx), drop = False)
    xmask = xmask.where(ymask.notnull())
    ymask = ymask.where(xmask.notnull())
    data = input_ds.where(xmask.notnull())

    # Transform input data from 2D to 1D and remove empty pixels.
    xmask_pixel = xmask.stack({"pixel": ("nx", "ny")}).dropna("pixel")
    ymask_pixel = ymask.stack({"pixel": ("nx", "ny")}).dropna("pixel")
    # NOTE: not using .dropna() because data needs to have same dimensions as `ymask_pixel`
    # and `xmask_pixel`, but can contain nan inside the chunk domain.
    data_pixel = data.stack({"pixel": ("nx", "ny")}).sel(pixel = xmask_pixel.pixel)

    # Load input data coordinates for scipy.
    xy = np.dstack((xmask_pixel.values,
                    ymask_pixel.values))[0]

    # Check if there is enough input data for current chunk.
    if xy.size < 20:
        return output_ds.unstack()

    # Determine distances between grid pixel and nearest input data point.
    tree = cKDTree(xy)
    xi = np.stack(np.meshgrid(grid_ds.x, grid_ds.y), axis = 2)
    dists = xr.DataArray(tree.query(xi)[0], dims = ("y", "x"), coords = { "y": grid_ds.y, "x": grid_ds.x})
    dists_pixel = dists.stack({"grid_pixel": ("y", "x")}).drop("grid_pixel")

    # Load pixel coordinates for scipy functions, excluding pixels 
    # for which input data is too far away.
    max_dist = np.min([max_px_dist * dx, max_px_dist * dy])
    grid_adj_ds = grid_adj_ds.where(dists_pixel < max_dist, drop = True)
    uv = np.dstack((grid_adj_ds.x.values, grid_adj_ds.y.values))[0]

    # 2D Delaunay tessellation.
    # Also see: https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids
    tri = Delaunay(xy)

    # Find simplex for grid points.
    simplex = tri.find_simplex(uv)

    # Determine vertices and weights.
    vtx = np.take(tri.simplices, simplex, axis = 0)
    temp = np.take(tri.transform, simplex, axis = 0)
    delta = uv - temp[:, 2]
    bary = np.einsum('njk,nk->nj', temp[:, :2, :], delta)
    wts = np.hstack((bary, 1 - bary.sum(axis = 1, keepdims = True)))

    # Wrap vertices and weights in xr.DataArray.
    wts_da = xr.DataArray(wts, dims  = ("grid_pixel", "j"), coords = {"j": range(3), "grid_pixel": grid_adj_ds.grid_pixel})
    vtx_da = xr.DataArray(vtx, dims  = ("grid_pixel", "j"), coords = {"j": range(3), "grid_pixel": grid_adj_ds.grid_pixel})

    # Select relevant input data.
    vls = data_pixel.isel(pixel = vtx_da)

    # Apply weights to input data.
    for var in input_ds.data_vars:
        da = xr.dot(vls[var], wts_da, dims = "j").where((wts_da > 0).all("j"))
        output_ds[var] = da.reindex(grid_pixel = dists_pixel.grid_pixel).drop("grid_pixel")

    return output_ds.unstack()

if __name__ == "__main__":

    # Small test.
    input_ds = xr.tutorial.open_dataset("rasm")
    input_ds = input_ds.rename_dims({"x": "nx", "y": "ny"}).rename_vars({"xc":"x", "yc": "y"})

    grid_ds = create_grid(input_ds, 3.0, 2.5).chunk({"x":500, "y":500})

    with ProgressBar():
        out = xr.map_blocks(regrid, grid_ds, (input_ds,), template = grid_ds).compute()

    from cartopy import crs as ccrs

    fig = plt.figure(figsize=(15, 10))
    ax1 = plt.subplot(2, 1, 1, projection = ccrs.PlateCarree())
    ax2 = plt.subplot(2, 1, 2, projection = ccrs.PlateCarree())
    kwargs = {"x": "x", "y": "y", "vmin": -20, "vmax": 20}
    bla = input_ds.isel(time=0).Tair.plot.pcolormesh(ax = ax1, **kwargs)
    ax1.coastlines()
    ax1.set_facecolor("lightgray")
    ax1.gridlines(draw_labels = True, color='gray', linestyle=':')
    ax1.set_title("curvilinear")
    out.isel(time=0).Tair.plot.pcolormesh(ax = ax2, **kwargs)
    ax2.coastlines()
    ax2.set_facecolor("lightgray")
    ax2.gridlines(draw_labels = True, color='gray', linestyle=':')
    ax2.set_title("rectolinear")
