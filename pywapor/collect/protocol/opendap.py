import os
import rasterio
import rioxarray.merge
import tempfile
import warnings
import requests
import copy
import urllib.parse
import xarray as xr
import numpy as np
from pywapor.general.logger import log
from rasterio.crs import CRS
from bs4 import BeautifulSoup
from urllib.parse import urlsplit, urlunsplit
from pywapor.enhancers.apply_enhancers import apply_enhancer
from pywapor.collect.protocol.crawler import download_url, download_urls
from pywapor.general.processing_functions import save_ds, process_ds, remove_ds

def download(fp, product_name, coords, variables, post_processors, 
                fn_func, url_func, un_pw = None, tiles = None,  
                data_source_crs = None, parallel = False, spatial_tiles = True, 
                request_dims = True, timedelta = None):
    """Download data from a OPENDaP server.

    Parameters
    ----------
    fp : str
        Path to file in which to download.
    product_name : str
        Name of product.
    coords : dict
        Coordinate names and boundaries.
    variables : dict
        Keys are variable names, values are additional settings.
    post_processors : dict
        Processors to apply to specific variables.
    url_func : function
        Function that takes `product_name` as input and return a url.
    un_pw : tuple, optional
        Username and password to use, by default None.
    tiles : list, optional
        Tiles to download, by default None.
    data_source_crs : rasterio.CRS.crs, optional
        CRS of datasource, by default None.
    parallel : bool, optional
        Download files in parallel (currently not implemented), by default False.
    spatial_tiles : bool, optional
        Whether the tiles are spatial or temporal, by default True.
    request_dims : bool, optional
        Include dimension settings in the OPENDaP request, by default True.
    timedelta : datetime.datetime.timedelta, optional
        Shift the time coordinates by tdelta, by default None.

    Returns
    -------
    xr.Dataset
        Dataset with downloaded data.
    """

    folder = os.path.split(fp)[0]

    # Create selection object.
    selection = create_selection(coords, target_crs = data_source_crs)

    # Make output filepaths, should be same length as `urls`.
    fps = [os.path.join(folder, fn_func(product_name, x)) for x in tiles]

    # Make data request URLs.
    session = start_session(url_func(product_name, tiles[0]), selection, un_pw)
    if spatial_tiles:
        idxss = [find_idxs(url_func(product_name, x), selection, session) for x in tiles]
        urls = [create_url(url_func(product_name, x), idxs, variables, request_dims = request_dims) for x, idxs in zip(tiles, idxss)]
    else:
        idxs = find_idxs(url_func(product_name, tiles[0]), selection, session)
        urls = [create_url(url_func(product_name, x), idxs, variables, request_dims = request_dims) for x in tiles]

    # Download data.
    files = download_urls(urls, "", session, fps = fps, parallel = parallel)

    # Merge spatial tiles.
    coords_ = {k: [v[0], selection[v[0]]] for k,v in coords.items()}
    if spatial_tiles:
        dss = [process_ds(xr.open_dataset(x, decode_coords = "all"), coords_, variables, crs = data_source_crs) for x in files]
        ds = rioxarray.merge.merge_datasets(dss)
    else:
        dss = files
        ds = xr.concat([xr.open_dataset(x, decode_coords="all") for x in files], dim = "time")
        ds = process_ds(ds, coords_, variables, crs = data_source_crs)

    # Reproject if necessary.
    if ds.rio.crs.to_epsg() != 4326:
        ds = ds.rio.reproject(CRS.from_epsg(4326))

    ds = ds.rio.clip_box(coords["x"][1][0], coords["y"][1][0], coords["x"][1][1], coords["y"][1][1])

    # Apply product specific functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    if isinstance(timedelta, np.timedelta64):
        ds["time"] = ds["time"] + timedelta

    # Remove unrequested variables.
    ds = ds[list(post_processors.keys())]
    
    # Save final output.
    ds.attrs = {}
    for var in ds.data_vars: # NOTE Keeping band attributes can cause problems when 
        # opening the data using rasterio (in reproject_chunk), see https://github.com/rasterio/rasterio/discussions/2751
        ds[var].attrs = {}

    ds = save_ds(ds, fp, encoding = "initiate", label = "Saving merged data.")

    # Remove temporary files.
    if not isinstance(dss, type(None)):
        for x in dss:
            remove_ds(x)

    return ds

def find_idxs(base_url, selection, session):
    fp = tempfile.NamedTemporaryFile(suffix=".nc").name
    url_coords = base_url + urllib.parse.quote(",".join(selection.keys()))
    fp = download_url(url_coords, fp, session, waitbar = False)
    ds = xr.open_dataset(fp, decode_coords = "all")
    idxs = dict()
    for k, v in selection.items():
        all_idxs = np.where((ds[k] >= v[0]) & (ds[k] <= v[1]))[0]
        if len(all_idxs) == 0:
            msg = f"--> No data found for `{k}` in range [{v[0]}, ...,{v[1]}]. Data is available between [{ds[k][0].values}, ..., {ds[k][-1].values}]."
            log.error(msg)
            raise ValueError(msg)
        idxs[k] = [np.min(all_idxs), np.max(all_idxs)]
    return idxs

def create_url(base_url, idxs, variables, request_dims = True):
    if request_dims:
        dims = [f"{k}[{v[0]}:{v[1]}]" for k, v in idxs.items()]
    else:
        dims = []
    varis = [f"{k}{''.join([f'[{idxs[dim][0]}:{idxs[dim][1]}]' for dim in v[0]])}" for k, v in variables.items()]
    url = base_url + urllib.parse.quote(",".join(dims + varis))
    return url

def start_session(base_url, selection, un_pw = [None, None]):
    if un_pw == [None, None]:
        warnings.filterwarnings("ignore", "password was not set. ")
    url_coords = base_url + urllib.parse.quote(",".join(selection.keys()))
    vrfy = {"NO": False, "YES": True}.get(os.environ.get("PYWAPOR_VERIFY_SSL", "YES"), True)
    session = setup_session('https://urs.earthdata.nasa.gov', username = un_pw[0], password = un_pw[1], check_url=url_coords, verify=vrfy)
    return session

def setup_session(uri,
                  username=None,
                  password=None,
                  check_url=None,
                  session=None,
                  verify=True,
                  username_field='username',
                  password_field='password'):
    '''
    A general function to set-up requests session with cookies
    using beautifulsoup and by calling the right url.
    '''

    if session is None:
        # Connections must be closed since some CAS
        # will cough when connections are kept alive:
        headers = [
            # ('User-agent', 'pydap/{}'.format(lib.__version__)),
                   ('Connection', 'close')]
        session = requests.Session()
        session.headers.update(headers)

    if uri is None:
        return session

    if not verify:
        verify_flag = session.verify
        session.verify = False

    if isinstance(uri, str):
        url = uri
    else:
        url = uri(check_url)

    if password is None or password == '':
        warnings.warn('password was not set. '
                      'this was likely unintentional '
                      'but will result is much fewer datasets.')
        if not verify:
            session.verify = verify_flag
        return session

    # Allow for several subsequent security layers:
    full_url = copy.copy(url)
    if isinstance(full_url, list):
        url = full_url[0]

    with warnings.catch_warnings():

        response = soup_login(session, url, username, password,
                              username_field=username_field,
                              password_field=password_field)

        # If there are further security levels.
        # At the moment only used for CEDA OPENID:
        if (isinstance(full_url, list) and
           len(full_url) > 1):
            for url in full_url[1:]:
                response = soup_login(session, response.url,
                                      username, password,
                                      username_field=None,
                                      password_field=None)
        response.close()

        if check_url:
            if (username is not None and
               password is not None):
                res = session.get(check_url, auth=(username, password))
                if res.status_code == 401:
                    res = session.get(res.url, auth=(username, password))
                res.close()
            raise_if_form_exists(check_url, session)

    if not verify:
        session.verify = verify_flag
    return session

def raise_if_form_exists(url, session):
    """
    This function raises a UserWarning if the link has forms
    """

    user_warning = ('Navigate to {0}, '.format(url) +
                    'login and follow instructions. '
                    'It is likely that you have to perform some one-time '
                    'registration steps before acessing this data.')

    resp = session.get(url)
    soup = BeautifulSoup(resp.content, 'lxml')
    if len(soup.select('form')) > 0:
        raise UserWarning(user_warning)


def soup_login(session, url, username, password,
               username_field='username',
               password_field='password'):
    resp = session.get(url)

    soup = BeautifulSoup(resp.content, 'lxml')
    login_form = soup.select('form')[0]

    def get_to_url(current_url, to_url):
        split_current = urlsplit(current_url)
        split_to = urlsplit(to_url)
        comb = [val2 if val1 == '' else val1
                for val1, val2 in zip(split_to, split_current)]
        return urlunsplit(comb)
    to_url = get_to_url(resp.url, login_form.get('action'))

    session.headers['Referer'] = resp.url

    payload = {}
    if username_field is not None:
        if len(login_form.findAll('input', {'name': username_field})) > 0:
            payload.update({username_field: username})

    if password_field is not None:
        if len(login_form.findAll('input', {'name': password_field})) > 0:
            payload.update({password_field: password})
        else:
            # If there is no password_field, it might be because
            # something should be handled in the browser
            # for the first attempt. This is common when using
            # pydap with the ESGF for the first time.
            raise Exception('Navigate to {0}. '
                            'If you are unable to '
                            'login, you must either '
                            'wait or use authentication '
                            'from another service.'
                            .format(url))

    # Replicate all other fields:
    for input in login_form.findAll('input'):
        if (input.get('name') not in payload and
           input.get('name') is not None):
            payload.update({input.get('name'): input.get('value')})

    # Remove other submit fields:
    submit_type = 'submit'
    submit_names = [input.get('name') for input
                    in login_form.findAll('input', {'type': submit_type})]
    for input in login_form.findAll('input', {'type': submit_type}):
        if ('submit' in submit_names and
           input.get('name').lower() != 'submit'):
            payload.pop(input.get('name'), None)

    return session.post(to_url, data=payload)

def download_xarray(url, fp, coords, variables, post_processors, 
                    data_source_crs = None, timedelta = None):
    """Download a OPENDaP dataset using xarray directly.

    Parameters
    ----------
    url : str
        URL to dataset.
    fp : str
        Path to file to download into.
    coords : dict
        Coordinates to request.
    variables : dict
        Variables to request.
    post_processors : dict
        Additional functions to apply to variables.
    data_source_crs : rasterio.CRS.crs, optional
        CRS of the data source, by default None.
    timedelta : datetime.datetime.timedelta, optional
        Shift the time coordinates by tdelta, by default None.

    Returns
    -------
    xr.Dataset
        Downloaded dataset.
    """

    warnings.filterwarnings("ignore", category=xr.SerializationWarning)
    online_ds = xr.open_dataset(url, decode_coords="all")
    # warnings.filterwarnings("default", category=xr.SerializationWarning)

    # Define selection.
    selection = create_selection(coords, target_crs = data_source_crs)

    # Make the selection on the remote.
    online_ds = online_ds.sel({k: slice(*v) for k, v in selection.items()})

    # Rename variables and assign crs.
    online_ds = process_ds(online_ds, coords, variables, crs = data_source_crs)

    # Download the data.
    ds = save_ds(online_ds, fp.replace(".nc", "_temp.nc"), label = f"Downloading data.")

    # Apply product specific functions.
    for var, funcs in post_processors.items():
        for func in funcs:
            ds, label = apply_enhancer(ds, var, func)
            log.info(label)

    if isinstance(timedelta, np.timedelta64):
        ds["time"] = ds["time"] + timedelta

    # Save final output
    out = save_ds(ds, fp, encoding = "initiate", label = "Saving netCDF.")

    remove_ds(ds)

    return out

def create_selection(coords, target_crs = None, source_crs = CRS.from_epsg(4326)):
    """Create a dictionary that can be given to `xr.Dataset.sel`.

    Parameters
    ----------
    coords : dict
        Dictionary describing the different dimensions over which to select. Possible keys
        are "x" for latitude, "y" for longitude and "t" for time, but other selections
        keys are also allowed (e.g. so select a band). Values are tuples with the first
        value the respective dimension names in the `ds` and the second value the selector.
    target_crs : rasterio.crs.CRS, optional
        crs of the dataset on which the selection will be applied, by default None.
    source_crs : rasterio.crs.CRS, optional
        crs of the `x` and `y` limits in `coords`, by default `epsg:4326`.

    Returns
    -------
    dict
        Dimension names with slices to apply to each dimension.
    """
    selection = {}

    if not isinstance(target_crs, type(None)):
        bounds = rasterio.warp.transform_bounds(source_crs, target_crs, 
                                                coords["x"][1][0], 
                                                coords["y"][1][0], 
                                                coords["x"][1][1], 
                                                coords["y"][1][1])
    else:
        bounds = [coords["x"][1][0], coords["y"][1][0], coords["x"][1][1], coords["y"][1][1]]
    
    selection[coords["x"][0]] = [bounds[0], bounds[2]]
    selection[coords["y"][0]] = [bounds[1], bounds[3]]

    if "t" in coords.keys():
        selection[coords["t"][0]] = [np.datetime64(t) for t in coords["t"][1]]

    for name, lim in coords.values():
        if name not in selection.keys():
            selection[name] = lim

    return selection
