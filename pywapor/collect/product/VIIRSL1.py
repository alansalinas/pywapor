#%%
import xarray as xr
from pywapor.general.processing_functions import save_ds
import os
from pywapor.general.curvilinear import create_grid, regrid
import datetime
import glob
workdir = r"/Users/hmcoerver/On My Mac/viirs_test"

latlim = [28.9, 29.7]
lonlim = [30.2, 31.2]

# Reformat bounding-box.
bb = [lonlim[0], latlim[0], lonlim[1], latlim[1]]

# Search for VNP02 images in workdir.
ncs02 = glob.glob(os.path.join(workdir, "VNP02IMG.A*.nc"))

# Make inventory of complete scenes in workdir.
scenes = dict()
for nc02 in ncs02:
    dt_str = ".".join(os.path.split(nc02)[-1].split(".")[1:3])
    ncs03 = glob.glob(os.path.join(workdir, f"VNP03IMG.{dt_str}*.nc"))
    if len(ncs03) != 1:
        continue
    else:
        nc03 = ncs03[0]
    dt = datetime.datetime.strptime(dt_str, "A%Y%j.%H%M")
    scenes[dt] = (nc02, nc03)

# Create list to store xr.Datasets from a single scene.
dss = list()

# Loop over the scenes.
for dt, (nc02, nc03) in scenes.items():

        # Open the datasets.
        ds1 = xr.open_dataset(nc02, group = "observation_data", decode_cf = False)
        ds2 = xr.open_dataset(nc03, group = "geolocation_data", decode_cf = False, chunks = "auto")

        # Convert DN to Brightness Temperature using the provided loopup table.
        bt_da = ds1.I05_brightness_temperature_lut.isel(number_of_LUT_values = ds1.I05)

        # Rename some things.
        ds = ds2[["latitude", "longitude"]].rename({"latitude": "y", "longitude": "x"})

        # Chunk and mask invalid pixels.
        ds["bt"] = bt_da.chunk("auto").where((bt_da >= bt_da.valid_min) & 
                                             (bt_da <= bt_da.valid_max) & 
                                             (bt_da != bt_da._FillValue))

        # Move `x` and `y` from variables to coordinates
        ds = ds.set_coords(["x", "y"])

        # Create mask for selecting the bounding-box in the data.
        buffer = 0.2
        mask = ((ds.y >= bb[1] - buffer) &
                (ds.y <= bb[3] + buffer) & 
                (ds.x >= bb[0] - buffer) & 
                (ds.x <= bb[2] + buffer))

        # Apply the mask.
        ds = ds.where(mask, drop = True)

        # Save intermediate file.
        fp = os.path.join(workdir, "temp.nc")
        ds = save_ds(ds, fp)

        # Create rectolinear grid.
        grid_ds = create_grid(ds, 0.0033, 0.0033)

        # Regrid from curvilinear to rectolinear grid.
        out = regrid(grid_ds, ds)

        # Set some metadata.
        out = out.rio.write_crs(4326)
        out = out.rio.clip_box(*bb)
        out.bt.attrs = {k: v for k, v in ds.bt.attrs.items() if k in ["long_name", "units"]}

        # Add time dimension.
        out = out.expand_dims({"time": 1}).assign_coords({"time": [dt]})

        dss.append(out)


# %%
