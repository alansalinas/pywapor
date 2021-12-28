import pywapor.collect as c

def collect_sources(param, sources, dl_args, extra_source_locations = None):
    files = list()
    for source in sources:

        if source == "PROBAV" and param == "NDVI":
            files.append(c.PROBAV.NDVI(**dl_args))
        elif source == "MOD13" and param == "NDVI":
            files.append(c.MOD13.NDVI(**dl_args))
        elif source == "MYD13" and param == "NDVI":
            files.append(c.MYD13.NDVI(**dl_args))

        elif source == "PROBAV" and param == "ALBEDO":
            files.append(c.PROBAV.ALBEDO(**dl_args))
        elif source == "MCD43" and param == "ALBEDO":
            files.append(c.MCD43.ALBEDO(**dl_args))

        elif source == "MOD11" and param == "LST":
            files.append(c.MOD11.LST(**dl_args))
        elif source == "MYD11" and param == "LST":
            files.append(c.MYD11.LST(**dl_args))

        elif source == "CHIRPS" and param == "PRECIPITATION":
            files.append(c.CHIRPS.PRECIPITATION(**dl_args))
        elif source == "SRTM" and param == "DEM":
            files.append([c.SRTM.DEM(**dl_args)])

        elif source == "GLOBCOVER" and param == "LULC":
            files.append([c.GLOBCOVER.LULC(**dl_args)])
        elif source == "WAPOR" and param == "LULC":
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
        elif source == "GEOS5" and param == "p_air_24_0":
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
        elif source == "MERRA2" and param == "p_air_24_0": # TODO should be "p_air_0_24"
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

        elif source == "MERRA2" and param == "SOLAR_RADIATION":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['swgnet']))

        else:
            sideload_files = c.sideloader.search_product_files(source, extra_source_locations[source])
            if len(sideload_files) > 0:
                files.append(sideload_files)

    return files

if __name__ == "__main__":

    import os
    import pywapor

    project_folder = r"/Users/hmcoerver/pywapor_notebooks"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    startdate = "2021-07-01"
    enddate = "2021-07-11"

    dl_args = {
            "Dir": os.path.join(project_folder, "RAW"), 
            "latlim": latlim, 
            "lonlim": lonlim, 
            "Startdate": startdate, 
            "Enddate": enddate,
            }

    raw_lst_files = collect_sources("LST", ['MOD11', 'MYD11'], dl_args, None)
    ds_lst = pywapor.pre_se_root.combine_lst(raw_lst_files)

    sds, eds, prds = pywapor.pre_se_root.calc_periods(ds_lst.time.values, "3H")

    dl_args["Startdate"] = sds
    dl_args["Enddate"] = eds
    dl_args["Periods"] = prds

    MERRAmeteo_vars = ['t2m', 'u2m', 'v2m', 'q2m',  'tpw', 'ps', 'slp']
    GEOSmeteo_vars = ['t2m', 'u2m', 'v2m', 'qv2m', 'tqv', 'ps', 'slp']

    # all_files = dict()
    for var in ["t_air_i", 
                "u2m_i", "v2m_i", "qv_i", "wv_i", "p_air_i", "p_air_0_i"
                ]:

        files = collect_sources(var, ["GEOS5", "MERRA2"], dl_args, extra_source_locations = None)

    # print(files)