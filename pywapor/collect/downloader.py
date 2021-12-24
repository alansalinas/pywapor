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
        elif source == "MERRA2" and param == "p_air_24_0":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['slp']))

        elif source == "MERRA2" and param == "SOLAR_RADIATION":
            files.append(c.MERRA2.daily_MERRA2(**dl_args, Vars = ['swgnet']))

        else:
            sideload_files = c.sideloader.search_product_files(source, extra_source_locations[source])
            if len(sideload_files) > 0:
                files.append(sideload_files)
    return files