import os
from dask.diagnostics import ProgressBar
import numpy as np
from pywapor.general.logger import log
import xarray as xr
from scipy.interpolate import griddata
import numpy as np
from scipy.spatial import cKDTree
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
import shutil
import glob
import rasterio.warp

def process_ds(ds, coords, variables, crs = None):

    ds = ds[list(variables.keys())]

    ds = ds.rename({v[0]:k for k,v in coords.items() if k in ["x", "y"]})
    ds = ds.rename({k: v[1] for k, v in variables.items()})

    if not isinstance(crs, type(None)):
        ds = ds.rio.write_crs(crs)

    ds = ds.rio.write_grid_mapping("spatial_ref")

    for var in [x for x in list(ds.variables) if x not in ds.coords]:
        if "grid_mapping" in ds[var].attrs.keys():
            del ds[var].attrs["grid_mapping"]

    ds = ds.sortby("y", ascending = False)
    ds = ds.sortby("x")

    ds.attrs = {}

    return ds

def save_ds(ds, fp, decode_coords = "all", encoding = None, chunks = "auto"):
    """Save a `xr.Dataset` as netcdf.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to save.
    fp : str
        Path to file to create.
    decode_coords : str, optional
        Controls which variables are set as coordinate variables when
        reopening the dataset, by default None.

    Returns
    -------
    xr.Dataset
        The newly created dataset.
    """
    temp_fp = fp.replace(".nc", "_temp")

    folder = os.path.split(fp)[0]
    if not os.path.isdir(folder):
        os.makedirs(folder)

    if isinstance(chunks, dict):
        chunks = {dim: v for dim, v in chunks.items() if dim in ds.dims}

    ds = ds.chunk(chunks)

    with ProgressBar(minimum = 50, dt = 2.0):
        ds.to_netcdf(temp_fp, engine = "netcdf4", encoding = encoding)

    ds = ds.close()

    os.rename(temp_fp, fp)

    ds = open_ds(fp, decode_coords = decode_coords, chunks = chunks)

    return ds

def open_ds(fp, decode_coords = "all", chunks = "auto"):
    ds = xr.open_dataset(fp, decode_coords = decode_coords, chunks = chunks)
    return ds

def create_wkt(latlim, lonlim):
    left = lonlim[0]
    bottom = latlim[0]
    right = lonlim[1]
    top = latlim[1]
    x = f"{left} {bottom},{right} {bottom},{right} {top},{right} {bottom},{left} {bottom}"
    return "GEOMETRYCOLLECTION(POLYGON((" + x + ")))"

def unpack(file, folder):
    fn = os.path.splitext(file)[0]
    shutil.unpack_archive(os.path.join(folder, file), folder)
    folder = [x for x in glob.glob(os.path.join(folder, fn + "*")) if os.path.isdir(x)][0]
    return folder

def transform_bb(src_crs, dst_crs, bb):
    bb =rasterio.warp.transform_bounds(src_crs, dst_crs, *bb, densify_pts=21)
    return bb

# def regrid_curvilinear(ds, resolution):

#     # Define new grid.
#     xmin = ds.x.min().values
#     xmax = ds.x.max().values
#     ymin = ds.y.min().values
#     ymax = ds.y.max().values
#     dx = dy = resolution
#     gridx = np.arange(xmin, xmax + dx, dx)
#     gridy = np.arange(ymin, ymax + dy, dy)
#     mgridx, mgridy = np.meshgrid(gridx, gridy)

#     out_ds = xr.Dataset(None, coords = {"y": gridy, "x": gridx})

#     for var in ds.data_vars:

#         values = ds[var].values.ravel()

#         # Remove no-data pixels.
#         y = ds.y.values.ravel()[np.isfinite(values)]
#         x = ds.x.values.ravel()[np.isfinite(values)]
#         values = values[np.isfinite(values)]
#         points = np.dstack((x,y))[0]

#         # Calculate pixel distances to points.
#         tree = cKDTree(points)
#         xi = _ndim_coords_from_arrays((mgridx, mgridy))
#         dists = tree.query(xi)[0]

#         # Interpolate points to grid.
#         grid_z0 = griddata(points, values, (mgridx, mgridy), method = "linear")
        
#         # Mask pixels too far away from any point.
#         grid_z0[dists > dx] = np.nan

#         # Wrap output into xr.Dataset.
#         out_ds[var] = xr.DataArray(grid_z0, coords = out_ds.coords)

#     out_ds = out_ds.rio.write_crs("epsg:4326")

#     out_ds = out_ds.sortby(out_ds.y, ascending=False)

#     return out_ds

# def ds_remove_except(ds, keep_vars):
#     """Remove all variables from a dataset except the variables specified with
#     `keep_vars`. Variables that are coordinates are never removed.

#     Parameters
#     ----------
#     ds : xr.Dataset
#         Dataset from which to remove variables.
#     keep_vars : list
#         List of variables to keep.

#     Returns
#     -------
#     xr.Dataset
#         Dataset from which variables have been removed.
#     """
#     keep_vars_coords = list()
#     for var in keep_vars:
#         if var in list(ds.variables):
#             keep_vars_coords += [k for k, v in ds[var].coords.items() if v.size > 0]
#     keep_vars_coords = keep_vars + np.unique(keep_vars_coords).tolist()
#     drop_vars = [x for x in list(ds.variables) if x not in keep_vars_coords]
#     ds = ds.drop_vars(drop_vars)
#     return ds

# def domain_overlaps_domain(domain1, domain2, partially = True):
#     """Chechs if `domain1` is (partically) contained inside `domain2`.

#     Parameters
#     ----------
#     domain1 : list
#         Borders of first domain.
#     domain2 : _type_
#         Borders of second domain.
#     partially : bool, optional
#         Whether or not `domain1` needs the be entirely inside `domain2`, by default True.

#     Returns
#     -------
#     bool
#         Whether `domain1` is (partially) overlaping with `domain2`.
#     """
#     check1 = bool(domain2[0] <= domain1[0] <= domain2[1])
#     check2 = bool(domain2[0] <= domain1[1] <= domain2[1])
#     if partially:
#         return check1 or check2
#     else:
#         return check1 and check2

# def reproject_ds(ds, fp, target_crs, source_crs = None):
#     # Remove unused coordinates.
#     ds_proj = ds.drop_vars([x for x in ds.coords if len(ds[x].dims) > 1])

#     # Remove existing grid_mapping data.
#     # for var in list(ds_proj.variables):
#     #     if "grid_mapping" in ds_proj[var].attrs.keys():
#     #         del ds_proj[var].attrs["grid_mapping"]

#     # Assign crs.
#     if not isinstance(source_crs, type(None)):
#         if isinstance(ds_proj.rio.crs, type(None)):
#             log.info("--> Setting crs.")
#         else:
#             log.warn("--> Overwriting crs.")
#         ds_proj = ds_proj.rio.write_crs(source_crs)

#     # Reproject to new crs.
#     ds_proj = ds_proj.rio.reproject(target_crs)

#     # Save output.
#     ds = save_ds(ds_proj, fp, decode_coords = "all")
#     return ds

# def export_ds_to_tif(ds, vars, base_folder):
#     """Export selected `vars` from a xr.Dataset (`ds`) into `base_folder` as geotiffs.

#     Parameters
#     ----------
#     ds : xr.Dataset
#         Dataset to export.
#     vars : list
#         Variables in `ds` to export.
#     base_folder : str
#         Path to folder in which to store geotiffs.

#     Returns
#     -------
#     dict
#         Keys are variable names, values are lists with paths to the generated
#         geotiffs.
#     """

#     all_files = dict()

#     if isinstance(base_folder, type(None)):
#         base_folder = os.path.split(ds.encoding["source"])[0]
    
#     for var in vars:

#         if var not in list(ds.variables):
#             continue

#         if "lat" not in list(ds[var].coords) or "lon" not in list(ds[var].coords):
#             continue

#         all_files[var] = list()

#         folder = os.path.join(base_folder, f"{var}_output")
        
#         if not os.path.exists(folder):
#             os.makedirs(folder)

#         if "source" in ds[var].attrs.keys():
#             source = ds[var].source
#         else:
#             source = "-"

#         if "unit" in ds[var].attrs.keys():
#             unit = ds[var].unit
#         else:
#             unit = "-"

#         if "time" in list(ds[var].coords):
#             iterator = ds.time
#             instantaneous = True
#         elif "epoch" in list(ds[var].coords):
#             iterator = ds.epoch
#             instantaneous = False
#         else:
#             log.warning(f"! --> not exporting {var}")

#         for t in iterator:

#             if instantaneous:
#                 tlength = "inst"
#                 date_str = pd.Timestamp(t.values).strftime("%Y.%m.%d.%H.%M")
#                 array = ds[var].sel(time = t).values
#             else:
#                 tdelta = ds.epoch_ends.sel(epoch = t) - ds.epoch_starts.sel(epoch = t)
#                 tlength = tdelta.values.astype('timedelta64[D]').astype(int)
#                 date_str = pd.Timestamp(ds.epoch_starts.sel(epoch = t).values).strftime("%Y.%m.%d")
#                 array = ds[var].sel(epoch = t).values

#             fn = f"{var.replace('_','-')}_{source}_{unit}_{tlength}_{date_str}.tif"
#             fh = os.path.join(folder, fn)

#             Save_as_tiff(fh, array, ds.geotransform, ds.projection)

#             all_files[var].append(fh)

#     return all_files

# def select_template(fhs):
#     """Given a list of lists with paths to geotiffs, determines which one of 
#     all the geotiffs has the highest resolution. Note that this is not based on 
#     pixel size, so passing geotiffs with very different domains can give unexpected
#     results.

#     Parameters
#     ----------
#     fhs : list
#         List of lists with paths to geotiffs.

#     Returns
#     -------
#     str
#         First file with the highest resolution.
#     xr.Dataset
#         Dataset with the required `lat` and `lon` values.
#     tuple
#         The geotransform of the selected geotiff.
#     resolution
#         Pixel size in units of the geotiffs projection.
#     """
#     # Flatten a list-of-lists into a  single list.
#     fhs = [val for sublist in fhs for val in sublist]

#     # Determine amount of pixels in each file.
#     sizes = [gdal.Open(fh).RasterXSize * gdal.Open(fh).RasterYSize for fh in fhs]

#     # Find index of file with most pixels.
#     idx = np.argmax(sizes)

#     # Create resample info.
#     example_fh = fhs[idx]
#     example_ds = xr.open_dataset(example_fh).isel(band = 0).drop_vars(["band", "spatial_ref"]).rename({"x": "lon", "y": "lat"})
#     example_geoinfo = get_geoinfo(example_fh)

#     # Calculate pixel size in meters.
#     resolution = get_resolution(example_fh)
#     log.info(f"--> Resampling resolution is ~{resolution:.0f} meter.")

#     return example_fh, example_ds, example_geoinfo, resolution

# def get_resolution(example_filepath):
#     """Get the average resolution of the given geotiff.

#     Parameters
#     ----------
#     example_filepath : str
#         Path to geotiff.

#     Returns
#     -------
#     float
#         Pixel size in units of the geotiffs projection.
#     """
#     ds = gdal.Open(example_filepath)
#     geo_ex = ds.GetGeoTransform()
#     xsize = ds.RasterXSize
#     ysize = ds.RasterYSize
#     dlat, dlon = calc_dlat_dlon(geo_ex, xsize, ysize)
#     dem_resolution = (np.nanmean(dlon) + np.nanmean(dlat))/2
#     return dem_resolution

# def reproject_clip(source_fp, dest_fp = None, bb = None,
#                     compress = True, dstSRS = "epsg:4326"):
#     """Function that calls gdal.Warp(), can be used to clip and reproject 
#     a geotiff in one go.

#     Parameters
#     ----------
#     source_fp : str
#         Path to file to be warped.
#     dest_fp : [type], optional
#         Warped output file, overwrites input when `None` given, by default None.
#     bb : tuple, optional
#         Clip the geotiff the a bounding-box.
#         First item is latlim, second item is lonlim, by default None.
#     compress : bool, optional
#         Compress the data using DEFLATE, by default True.
#     dstSRS : str, optional
#         Reproject the output to `dstSRS`, by default "epsg:4326".
#     """
#     options_dict = {"dstSRS": dstSRS}
    
#     if not isinstance(bb, type(None)):
#         options_dict["outputBounds"] = (bb[1][0], bb[0][0], 
#                                         bb[1][1], bb[0][1])

#     if compress:
#         options_dict["creationOptions"] = ["COMPRESS=DEFLATE", "ZLEVEL=8"]

#     if isinstance(dest_fp, type(None)):
#         log.debug(f"Overwriting input file ({source_fp}).")
#         dest_fp = source_fp

#     options = gdal.WarpOptions(**options_dict)

#     out = gdal.Warp(dest_fp, source_fp, options = options)
#     out.FlushCache()
#     out = None

#     return dest_fp

# def apply_mask(a, indices, axis):
#     """Select indices in `a` along an axis.

#     Parameters
#     ----------
#     a : np.ndarray
#         Array with one dimension more than `indices`.
#     indices : np.ndarray
#         Array used to index `a`, should have one dimension fewer than `a`.
#     axis : int
#         Axis along which to select values from `a`.

#     Returns
#     -------
#     np.ndarray
#         Array with same shape as `indices`.
#     """
#     # https://stackoverflow.com/questions/15469302/numpy-3d-to-2d-transformation-based-on-2d-mask-array
#     # TODO replace with np.take(....)?

#     magic_index = [np.arange(i) for i in indices.shape]
#     magic_index = np.ix_(*magic_index)
#     magic_index = magic_index[:axis] + (indices,) + magic_index[axis:]
#     return a[magic_index]

# def reproj_file(file, template, method):
#     """Reprojects and resamples a file and return the output as an np.ndarray.

#     Parameters
#     ----------
#     file : {str | gdal.Dataset}
#         File to reproject.
#     template : {str | gdal.Dataset}
#         File to use a reprojection template.
#     method : int
#         Interpolation method, see 'pywapor.processing_functions.reproject_dataset_example'
#         for more info.

#     Returns
#     -------
#     np.ndarray
#         The resamples data.
#     """
#     ds = reproject_dataset_example(file, template, method = method)
#     array = open_as_array(ds)
#     return array

# def combine_dicts(dicts):
#     """Combines dictionaries by appending values for identical keys in a list.

#     Parameters
#     ----------
#     dicts : list
#         A list of dictionaries.

#     Returns
#     -------
#     dict
#         Dictionary with all the values and keys from the dictionaries passed inside 
#         the list `dicts`.
#     """
#     new_dict = dict()
#     for d in dicts:
#         for key, value in d.items():
#             if key in new_dict.keys():
#                 new_dict[key].append(value)
#             else:
#                 new_dict[key] = [value]
#     return new_dict

# def get_geoinfo(template_file):
#     """Extract relevant geo information from a geotiff file.

#     Parameters
#     ----------
#     template_file : str
#         Path to geotiff file.

#     Returns
#     -------
#     tuple
#         Tuple with the geotransform, projection, amount of pixels in x-direction and
#         amount of pixels in y-direction.
#     """
#     if isinstance(template_file, str):
#         ds = gdal.Open(template_file)
#     elif isinstance(template_file, gdal.Dataset):
#         ds = template_file
#     geo_ex = ds.GetGeoTransform()
#     proj_ex = ds.GetProjection()
#     size_x_ex = ds.RasterXSize
#     size_y_ex = ds.RasterYSize
#     return (geo_ex, proj_ex, size_x_ex, size_y_ex)

# def Extract_Data_gz(zip_filename, outfilename):
#     """Extract a zip-file and removes the zip-file afterwards.

#     Parameters
#     ----------
#     zip_filename : str
#         Path to zip-file.
#     outfilename : str
#         Path to output.
#     """
#     # TODO maybe move to CHIRPS?
#     with gzip.GzipFile(zip_filename, 'rb') as zf:
#         file_content = zf.read()
#         save_file_content = open(outfilename, 'wb')
#         save_file_content.write(file_content)
#     save_file_content.close()
#     zf.close()
#     os.remove(zip_filename)
    
# def Save_as_MEM(data='', geo='', projection=''):
#     """Create a gdal.Dataset in memory from a np.ndarray.

#     Parameters
#     ----------
#     data : np.ndarray
#         Data to put inside the gdal.Dataset.
#     geo : list
#         Geotransform
#     projection : osr.SpatialReference
#         Projection to be used in the gdal.Dataset, by default "WGS84".
#     """
#     # TODO fix default values '' etc.
#     driver = gdal.GetDriverByName("MEM")
#     dst_ds = driver.Create('', int(data.shape[1]), int(data.shape[0]), 1,
#                            gdal.GDT_Float32)
#     srse = osr.SpatialReference()
#     if projection == '':
#         srse.SetWellKnownGeogCS("WGS84")

#     else:
#         # TODO Remove try/excepts
#         try:
#             if not srse.SetWellKnownGeogCS(projection) == 6:
#                 srse.SetWellKnownGeogCS(projection)
#             else:
#                 try:
#                     srse.ImportFromEPSG(int(projection))
#                 except:
#                     srse.ImportFromWkt(projection)
#         except:
#             try:
#                 srse.ImportFromEPSG(int(projection))
#             except:
#                 srse.ImportFromWkt(projection)
#     dst_ds.SetProjection(srse.ExportToWkt())
#     dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
#     dst_ds.SetGeoTransform(geo)
#     dst_ds.GetRasterBand(1).WriteArray(data)
#     return(dst_ds)   
    
# def Save_as_tiff(name='', data='', geo='', projection=''):
#     """Create a geotiff from a np.ndarray.

#     Parameters
#     ----------
#     nme : str
#         Path to output file.
#     data : np.ndarray
#         Data to put inside the geotiff.
#     geo : list
#         Geotransform
#     projection : osr.SpatialReference
#         Projection to be used in the geotiff, by default "WGS84".
#     """
#     # TODO fix default values '' etc.
#     # Change no data values
#     data[np.isnan(data)] = -9999

#     if not os.path.exists(os.path.split(name)[0]):
#         os.makedirs(os.path.split(name)[0])
    
#     # save as a geotiff
#     driver = gdal.GetDriverByName("GTiff")
#     dst_ds = driver.Create(name, int(data.shape[1]), int(data.shape[0]), 1,
#                            gdal.GDT_Float32, ['COMPRESS=LZW', 'PREDICTOR=3'])
#     srse = osr.SpatialReference()
#     if projection == '':
#         srse.SetWellKnownGeogCS("WGS84")

#     # Set the projection, which can be an EPSG code or a well known GeogCS
#     else:
#         # TODO get rid of all these try/excepts
#         try:
#             if not srse.SetWellKnownGeogCS(projection) == 6:
#                 srse.SetWellKnownGeogCS(projection)
#             else:
#                 try:
#                     srse.ImportFromEPSG(int(projection))
#                 except:
#                     srse.ImportFromWkt(projection)
#         except:
#             try:
#                 srse.ImportFromEPSG(int(projection))
#             except:
#                 srse.ImportFromWkt(projection)
                
#     # Save the tiff file
#     dst_ds.SetProjection(srse.ExportToWkt())
#     dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
#     dst_ds.SetGeoTransform(geo)
#     dst_ds.GetRasterBand(1).WriteArray(data)
#     dst_ds = None
    
#     return()
    
# def reproject_MODIS(input_name, epsg_to):
#     # TODO move this to pywapor.collect.MODIS
    
#     # Define the output name
#     name_out = ''.join(input_name.split(".")[:-1]) + '_reprojected.tif'
   
#     src_ds = gdal.Open(input_name)
    
#     # Define target SRS
#     dst_srs = osr.SpatialReference()
#     dst_srs.ImportFromEPSG(int(epsg_to))
#     dst_wkt = dst_srs.ExportToWkt()
    
#     error_threshold = 0.125  # error threshold --> use same value as in gdalwarp
#     resampling = gdal.GRA_NearestNeighbour
    
#     # Call AutoCreateWarpedVRT() to fetch default values for target raster dimensions and geotransform
#     tmp_ds = gdal.AutoCreateWarpedVRT( src_ds,
#                                    None, # src_wkt : left to default value --> will use the one from source
#                                    dst_wkt,
#                                    resampling,
#                                    error_threshold )
#     dst_ds = gdal.GetDriverByName('GTiff').CreateCopy(name_out, tmp_ds)
#     dst_ds = None 

#     return(name_out)
    
# def clip_data(input_file, latlim, lonlim):
#     """Clip the extents of a geotiff.

#     Parameters
#     ----------
#     input_file : {str | gdal.Dataset}
#         Dataset to be clipped
#     latlim : list
#         List with lower and upper latitude boundaries.
#     lonlim : list
#         List with lower and upper longitude boundaries.

#     Returns
#     -------
#     np.ndarray
#         The clipped data.
#     tuple
#         The new geotransform.
#     """
#     # TODO merge with reproject_clip
#     # TODO remove try/except.
#     try:
#         if input_file.split('.')[-1] == 'tif':
#             dest_in = gdal.Open(input_file)
#         else:
#             dest_in = input_file
#     except:
#         dest_in = input_file

#     # Open Array
#     data_in = dest_in.GetRasterBand(1).ReadAsArray()

#     # Define the array that must remain
#     Geo_in = dest_in.GetGeoTransform()
#     Geo_in = list(Geo_in)
#     Start_x = np.max([int(np.floor(((lonlim[0]) - Geo_in[0])/ Geo_in[1])),0])
#     End_x = np.min([int(np.ceil(((lonlim[1]) - Geo_in[0])/ Geo_in[1])),int(dest_in.RasterXSize)])

#     Start_y = np.max([int(np.floor((Geo_in[3] - latlim[1])/ -Geo_in[5])),0])
#     End_y = np.min([int(np.ceil(((latlim[0]) - Geo_in[3])/Geo_in[5])), int(dest_in.RasterYSize)])

#     #Create new GeoTransform
#     Geo_in[0] = Geo_in[0] + Start_x * Geo_in[1]
#     Geo_in[3] = Geo_in[3] + Start_y * Geo_in[5]
#     Geo_out = tuple(Geo_in)

#     data = np.zeros([End_y - Start_y, End_x - Start_x])

#     data = data_in[Start_y:End_y,Start_x:End_x]
#     dest_in = None

#     return(data, Geo_out)    
    
    
# def Extract_Data(input_file, output_folder):
#     """Extract a zipfile.

#     Parameters
#     ----------
#     input_file : str
#         Path to zip-file.
#     output_folder : str
#         Path to output folder.
#     """
#     # TODO merge with pywapor.general.processing_functions.Extract_Data_gz
#     # extract the data, used in SRTM and Globcover
#     z = zipfile.ZipFile(input_file, 'r')
#     z.extractall(output_folder)
#     z.close()    

# def reproj_ds(source_file, example_file, ndv = -9999, override_dtype = True):
#     """Reproject and resample (nearest neighbour) a geotiff to match with
#     `example_file`.

#     Parameters
#     ----------
#     source_file : str
#         Path to source file.
#     example_file : str
#         Path to example file.
#     ndv : int, optional
#         No-data-value to be used in the reprojected dataset, by default -9999
#     override_dtype : bool, optional
#         Use the same dtype as `example_file` (True) or as `source_file` (False), 
#         by default True.

#     Returns
#     -------
#     gdal.Dataset
#         The reprojected dataset.
#     """
#     # TODO merge with reproject_dataset_example
#     source_ds = gdal.Open(source_file)

#     example_ds = gdal.Open(example_file)
#     example_xsize = example_ds.RasterXSize
#     example_ysize = example_ds.RasterYSize

#     if override_dtype:
#         example_dtype = example_ds.GetRasterBand(1).DataType
#     else:
#         example_dtype = source_ds.GetRasterBand(1).DataType

#     mem_drv = gdal.GetDriverByName('MEM')
#     dest_ds = mem_drv.Create('', example_xsize, example_ysize, 1, example_dtype)
#     dest_ds.SetGeoTransform(example_ds.GetGeoTransform())
#     dest_ds.SetProjection(example_ds.GetProjection())
#     dest_ds.GetRasterBand(1).Fill(ndv)
#     dest_ds.GetRasterBand(1).SetNoDataValue(ndv)

#     gdal.ReprojectImage(source_ds,
#                         dest_ds,
#                         source_ds.GetProjection(),
#                         dest_ds.GetProjection(),
#                         gdal.GRA_NearestNeighbour)

#     return dest_ds

# def reproject_dataset_example(dataset, dataset_example, method=1):
#     """Reproject and resample a dataset to match with `dataset_example`.

#     Parameters
#     ----------
#     dataset : {str | gdal.Dataset}
#         The dataset to be reprojected.
#     dataset_example : {str | gdal.Dataset}
#         The dataset to be used as a reprojection example.
#     method : int, optional
#         Select which method to use for resampling, with:
#         1 = gdal.GRA_NearestNeighbour,
#         2 = gdal.GRA_Bilinear,
#         3 = gdal.GRA_Lanczos,
#         4 = gdal.GRA_Average,
#         5 = gdal.GRA_Cubic,
#         6 = gdal.GRA_CubicSpline,
#         7 = gdal.GRA_Mode,
#         8 = gdal.GRA_Max,
#         9 = gdal.GRA_Min,
#         10 = gdal.GRA_Med,
#         11 = gdal.GRA_Q1,
#         12 = gdal.GRA_Q3,
#         13 = gdal.GRA_Sum, by default 1.

#     Returns
#     -------
#     gdal.Dataset
#         The reprojected dataset.

#     """
#     # TODO merge with reproj_ds
#     # open dataset that must be transformed
#     if isinstance(dataset, str):
#         g = gdal.Open(dataset)
#     elif isinstance(dataset, gdal.Dataset):
#         g = dataset
#     else:
#         raise ValueError
#     epsg_from = Get_epsg(g)

#     #exceptions
#     if epsg_from == 9001:
#         epsg_from = 5070

#     # open dataset that is used for transforming the dataset
#     if isinstance(dataset_example, str):
#         gland = gdal.Open(dataset_example)
#     elif isinstance(dataset_example, gdal.Dataset):
#         gland = dataset_example
#     else:
#         raise ValueError
#     epsg_to = Get_epsg(gland)

#     # Set the EPSG codes
#     osng = osr.SpatialReference()
#     osng.ImportFromEPSG(epsg_to)
#     wgs84 = osr.SpatialReference()
#     wgs84.ImportFromEPSG(epsg_from)

#     # Get shape and geo transform from example
#     geo_land = gland.GetGeoTransform()
#     col=gland.RasterXSize
#     rows=gland.RasterYSize

#     # Create new raster
#     mem_drv = gdal.GetDriverByName('MEM')
#     dest = mem_drv.Create('', col, rows, 1, gdal.GDT_Float32)
#     dest.SetGeoTransform(geo_land)
#     dest.SetProjection(osng.ExportToWkt())

#     # https://gis.stackexchange.com/questions/158503/9999-no-data-value-becomes-0-when-writing-array-to-gdal-memory-file
#     ndv = g.GetRasterBand(1).GetNoDataValue()
#     band = dest.GetRasterBand(1)
#     band.SetNoDataValue(ndv)
#     band.Fill(ndv, 0.0)

#     # Perform the projection/resampling
#     if method == 1:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_NearestNeighbour)
#     if method == 2:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Bilinear)
#     if method == 3:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Lanczos)
#     if method == 4:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Average)
#     if method == 5:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Cubic)
#     if method == 6:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_CubicSpline)
#     if method == 7:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Mode)
#     if method == 8:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Max)
#     if method == 9:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Min)
#     if method == 10:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Med)
#     if method == 11:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Q1)
#     if method == 12:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Q3)
#     if method == 13:
#         gdal.ReprojectImage(g, dest, wgs84.ExportToWkt(), osng.ExportToWkt(), gdal.GRA_Sum)
        
#     return(dest)

# def Get_epsg(g):
#     # TODO deprecate
#     try:
#         # Get info of the dataset that is used for transforming
#         g_proj = g.GetProjection()
#         Projection=g_proj.split('EPSG","')
#         epsg_to=int((str(Projection[-1]).split(']')[0])[0:-1])            
#     except:
#         epsg_to=4326

#     return(epsg_to)

# def open_as_array(input):
#     """Open dataset as array, replacing no-data-values with np.nan.

#     Parameters
#     ----------
#     input : {str | gdal.Dataset}
#         The dataset to open.

#     Returns
#     -------
#     np.ndarray
#         The data contained inside `input`.
#     """
#     if isinstance(input, str):
#         ds = gdal.Open(input)
#     elif isinstance(input, gdal.Dataset):
#         ds = input
#     array = ds.GetRasterBand(1).ReadAsArray().astype(float)
#     ndv = ds.GetRasterBand(1).GetNoDataValue()
#     array[array == ndv] = np.nan
#     return array  

# def Open_tiff_array(filename='', band=''):
#     # TODO deprecate
#     f = gdal.Open(filename)
#     if f is None:
#         print('%s does not exists' %filename)
#     else:
#         if band == '':
#             band = 1
#         Data = f.GetRasterBand(band).ReadAsArray()
#     return(Data)

# def Open_array_info(filename=''):
#     # TODO deprecate
#     try:
#         if filename.split('.')[-1] == 'tif':
#             f = gdal.Open(r"%s" %filename)
#         else:
#             f = filename
#     except:
#             f = filename       
#     try:
#         geo_out = f.GetGeoTransform()
#         proj = f.GetProjection()
#         size_X = f.RasterXSize
#         size_Y = f.RasterYSize
#         f = None
#     except:
#         print('%s does not exists' %filename)
        
#     return(geo_out, proj, size_X, size_Y)

# def gap_filling(dataset, NoDataValue, method = 1):
#     # TODO deprecate
#     try:
#         if dataset.split('.')[-1] == 'tif':
#             # Open the numpy array
#             data = Open_tiff_array(dataset)
#             Save_as_tiff_bool = 1
#         else:
#             data = dataset
#             Save_as_tiff_bool = 0
#     except:
#         data = dataset
#         Save_as_tiff_bool = 0

#     # fill the no data values
#     if NoDataValue == np.nan:
#         mask = ~(np.isnan(data))
#     else:
#         mask = ~(data==NoDataValue)
#     xx, yy = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
#     xym = np.vstack( (np.ravel(xx[mask]), np.ravel(yy[mask])) ).T
#     data0 = np.ravel( data[:,:][mask] )

#     if method == 1:
#         interp0 = NearestNDInterpolator( xym, data0 )
#         data_end = interp0(np.ravel(xx), np.ravel(yy)).reshape( xx.shape )

#     if method == 2:
#         interp0 = LinearNDInterpolator( xym, data0 )
#         data_end = interp0(np.ravel(xx), np.ravel(yy)).reshape( xx.shape )

#     if Save_as_tiff_bool == 1:
#         EndProduct=dataset[:-4] + '_GF.tif'

#         # collect the geoinformation
#         geo_out, proj, size_X, size_Y = Open_array_info(dataset)

#         # Save the filled array as geotiff
#         Save_as_tiff(name=EndProduct, data=data_end, geo=geo_out, projection=proj)

#     else:
#         EndProduct = data_end

#     return (EndProduct)  
    
def calc_dlat_dlon(geo_out, size_X, size_Y, lat_lon = None):
    """
    Calculated the dimensions of each pixel in meter.

    Parameters
    ----------
    geo_out: list
        Geotransform function of the array.
    size_X: int
        Number of pixels in x-direction.
    size_Y: int
        Number of pixels in y-direction.

    Returns
    -------
    np.ndarray
        Size of every pixel in the y-direction in meters.
    dlon: array
        Size of every pixel in the x-direction in meters.
    """
    if isinstance(lat_lon, type(None)):
        # Create the lat/lon rasters
        lon = np.arange(size_X + 1)*geo_out[1]+geo_out[0] - 0.5 * geo_out[1]
        lat = np.arange(size_Y + 1)*geo_out[5]+geo_out[3] - 0.5 * geo_out[5]
    else:
        lat, lon = lat_lon

    dlat_2d = np.array([lat,]*int(np.size(lon,0))).transpose()
    dlon_2d =  np.array([lon,]*int(np.size(lat,0)))

    # Radius of the earth in meters
    R_earth = 6371000

    # Calculate the lat and lon in radians
    lonRad = dlon_2d * np.pi/180
    latRad = dlat_2d * np.pi/180

    # Calculate the difference in lat and lon
    lonRad_dif = abs(lonRad[:,1:] - lonRad[:,:-1])
    latRad_dif = abs(latRad[:-1] - latRad[1:])

    # Calculate the distance between the upper and lower pixel edge
    a = np.sin(latRad_dif[:,:-1]/2) * np.sin(latRad_dif[:,:-1]/2)
    clat = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    dlat = R_earth * clat

    # Calculate the distance between the eastern and western pixel edge
    b = np.cos(latRad[1:,:-1]) * np.cos(latRad[:-1,:-1]) * np.sin(lonRad_dif[:-1,:]/2) * np.sin(lonRad_dif[:-1,:]/2)
    clon = 2 * np.arctan2(np.sqrt(b), np.sqrt(1-b))
    dlon = R_earth * clon

    return(dlat, dlon)