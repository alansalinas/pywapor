import pywapor.collect as c
import tqdm
import requests
import os


def collect_sources(param, sources, dl_args, extra_source_locations = None):
    files = list()
    for source in sources:

        if source == "PROBAV" and param == "ndvi":
            files.append(c.PROBAV.NDVI(**dl_args))
        elif source == "MOD13" and param == "ndvi":
            files.append(c.MOD13.NDVI(**dl_args))
        elif source == "MYD13" and param == "ndvi":
            files.append(c.MYD13.NDVI(**dl_args))

        elif source == "PROBAV" and param == "r0":
            files.append(c.PROBAV.ALBEDO(**dl_args))
        elif source == "MCD43" and param == "r0":
            files.append(c.MCD43.ALBEDO(**dl_args))

        elif source == "MOD11" and param == "lst":
            files.append(c.MOD11.LST(**dl_args))
        elif source == "MYD11" and param == "lst":
            files.append(c.MYD11.LST(**dl_args))

        elif source == "CHIRPS" and param == "p_24":
            files.append(c.CHIRPS.PRECIPITATION(**dl_args))
        elif source == "SRTM" and param == "z":
            files.append([c.SRTM.DEM(**dl_args)])

        elif source == "GLOBCOVER" and param == "lulc":
            files.append([c.Globcover.LULC(**dl_args)])
        elif source == "WAPOR" and param == "lulc":
            files.append(c.WAPOR.LULC(**dl_args))

        elif source == "GEOS5" and param == "t_air_24":
            files.append(c.GEOS5.daily(**dl_args, Vars = ['t2m']))
        elif source == "GEOS5" and param == "t_air_max_24":
            files.append(c.GEOS5.daily(**dl_args, Vars = ['t2m-max']))
        elif source == "GEOS5" and param == "t_air_min_24":
            files.append(c.GEOS5.daily(**dl_args, Vars = ['t2m-min']))
        elif source == "GEOS5" and param == "u2m_24":
            files.append(c.GEOS5.daily(**dl_args, Vars = ['u2m']))
        elif source == "GEOS5" and param == "v2m_24":
            files.append(c.GEOS5.daily(**dl_args, Vars = ['v2m']))
        elif source == "GEOS5" and param == "qv_24":
            files.append(c.GEOS5.daily(**dl_args, Vars = ['qv2m']))
        elif source == "GEOS5" and param == "p_air_0_24":
            files.append(c.GEOS5.daily(**dl_args, Vars = ['slp']))

        elif source == "GEOS5" and param == "t_air_i":
            files.append(c.GEOS5.three_hourly(**dl_args, Vars = ["t2m"]))
        elif source == "GEOS5" and param == "u2m_i":
            files.append(c.GEOS5.three_hourly(**dl_args, Vars = ["u2m"]))
        elif source == "GEOS5" and param == "v2m_i":
            files.append(c.GEOS5.three_hourly(**dl_args, Vars = ["v2m"]))
        elif source == "GEOS5" and param == "qv_i":
            files.append(c.GEOS5.three_hourly(**dl_args, Vars = ["qv2m"]))
        elif source == "GEOS5" and param == "wv_i":
            files.append(c.GEOS5.three_hourly(**dl_args, Vars = ["tqv"]))
        elif source == "GEOS5" and param == "p_air_i":
            files.append(c.GEOS5.three_hourly(**dl_args, Vars = ["ps"]))
        elif source == "GEOS5" and param == "p_air_0_i":
            files.append(c.GEOS5.three_hourly(**dl_args, Vars = ["slp"]))

        elif source == "MERRA2" and param == "t_air_24":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['t2m']))
        elif source == "MERRA2" and param == "t_air_max_24":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['t2m'], data_type = ["max"]))
        elif source == "MERRA2" and param == "t_air_min_24":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['t2m'], data_type = ["min"]))
        elif source == "MERRA2" and param == "u2m_24":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['u2m']))
        elif source == "MERRA2" and param == "v2m_24":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['v2m']))
        elif source == "MERRA2" and param == "qv_24":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['q2m']))
        elif source == "MERRA2" and param == "p_air_0_24":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['slp']))

        elif source == "MERRA2" and param == "t_air_i":
            files.append(c.MERRA2.hourly_MERRA2(**dl_args, Vars = ["t2m"]))
        elif source == "MERRA2" and param == "u2m_i":
            files.append(c.MERRA2.hourly_MERRA2(**dl_args, Vars = ["u2m"]))
        elif source == "MERRA2" and param == "v2m_i":
            files.append(c.MERRA2.hourly_MERRA2(**dl_args, Vars = ["v2m"]))
        elif source == "MERRA2" and param == "qv_i":
            files.append(c.MERRA2.hourly_MERRA2(**dl_args, Vars = ["q2m"]))
        elif source == "MERRA2" and param == "wv_i":
            files.append(c.MERRA2.hourly_MERRA2(**dl_args, Vars = ["tpw"]))
        elif source == "MERRA2" and param == "p_air_i":
            files.append(c.MERRA2.hourly_MERRA2(**dl_args, Vars = ["ps"]))
        elif source == "MERRA2" and param == "p_air_0_i":
            files.append(c.MERRA2.hourly_MERRA2(**dl_args, Vars = ["slp"]))

        elif source == "MERRA2" and param == "ra_24":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['swgnet']))

        elif source == "STATICS" and param == "land_mask":
            files.append(c.STATICS.collect(**dl_args, vars = ['land_mask']))
        elif source == "STATICS" and param == "lw_offset":
            files.append(c.STATICS.collect(**dl_args, vars = ['lw_offset']))
        elif source == "STATICS" and param == "lw_slope":
            files.append(c.STATICS.collect(**dl_args, vars = ['lw_slope']))
        elif source == "STATICS" and param == "r0_bare":
            files.append(c.STATICS.collect(**dl_args, vars = ['r0_bare']))
        elif source == "STATICS" and param == "r0_full":
            files.append(c.STATICS.collect(**dl_args, vars = ['r0_full']))
        elif source == "STATICS" and param == "rn_offset":
            files.append(c.STATICS.collect(**dl_args, vars = ['rn_offset']))
        elif source == "STATICS" and param == "rn_slope":
            files.append(c.STATICS.collect(**dl_args, vars = ['rn_slope']))
        elif source == "STATICS" and param == "rs_min":
            files.append(c.STATICS.collect(**dl_args, vars = ['rs_min']))
        elif source == "STATICS" and param == "t_amp_year":
            files.append(c.STATICS.collect(**dl_args, vars = ['t_amp_year']))
        elif source == "STATICS" and param == "t_opt":
            files.append(c.STATICS.collect(**dl_args, vars = ['t_opt']))
        elif source == "STATICS" and param == "vpd_slope":
            files.append(c.STATICS.collect(**dl_args, vars = ['vpd_slope']))
        elif source == "STATICS" and param == "z_obst_max":
            files.append(c.STATICS.collect(**dl_args, vars = ['z_obst_max']))
        elif source == "STATICS" and param == "z_oro":
            files.append(c.STATICS.collect(**dl_args, vars = ['z_oro']))

        elif (source, param) in extra_source_locations.keys():
            sideload_files = c.sideloader.search_product_files(source, extra_source_locations[(source, param)])
            if len(sideload_files) > 0:
                files.append(sideload_files)
        else:
            raise ValueError

        # Remove sources for which no files were downloaded.
        files = [x for x in files if len(x) > 0]

    return files

def url_to_file(url, out_file):

    file_object = requests.get(url)
    file_object.raise_for_status()

    total_size = int(file_object.headers.get('content-length', 0))

    waitbar = tqdm.tqdm(total = total_size, unit='Bytes', unit_scale=True, position = 0)

    folder = os.path.split(out_file)[0]
    if not os.path.exists(folder):
        os.makedirs(folder)

    with open(out_file, 'wb') as z:
        for data in file_object.iter_content(chunk_size=1024):
            size = z.write(data)
            waitbar.update(size)

if __name__ == "__main__":

    project_folder = r"/Users/hmcoerver/pywapor_notebooks"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    startdate = "2021-07-01"
    enddate = "2021-07-03"

    dl_args = {
            "Dir": os.path.join(project_folder, "RAW"),
            "latlim": latlim,
            "lonlim": lonlim,
            "Startdate": startdate,
            "Enddate": enddate,
            }

    raw_r0_files = collect_sources("r0", ['MCD43'], dl_args, None)
    raw_lst_files = collect_sources("lst", ['MOD11', 'MYD11'], dl_args, None)
    raw_ndvi_files = collect_sources("ndvi", ['MYD13', 'MOD13'], dl_args, None)

    # ds_lst = pywapor.pre_se_root.combine_lst(raw_lst_files)

    # sds, eds, prds = pywapor.pre_se_root.calc_periods(ds_lst.time.values, "3H")

    # dl_args["Startdate"] = sds
    # dl_args["Enddate"] = eds
    # dl_args["Periods"] = prds

    # MERRAmeteo_vars = ['t2m', 'u2m', 'v2m', 'q2m',  'tpw', 'ps', 'slp']
    # GEOSmeteo_vars = ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp']

    # # all_files = dict()
    # for var in ["t_air_i", 
    #             "u2m_i", "v2m_i", "qv_i", "wv_i", "p_air_i", "p_air_0_i"
    #             ]:

    #     files = collect_sources(var, ["GEOS5", "MERRA2"], dl_args, extra_source_locations = None)

    # print(files)