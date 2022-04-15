
def collect_sources(folder, sources, latlim, lonlim, timelim, return_fps = True):
    
    reversed_sources = reverse_sources(sources)

    dss = dict()

    for (source, product_name), req_vars in reversed_sources.items():

        dl_module = __import__(f"pywapor.collect_new.product.{source}", 
                            fromlist=[source])

        args = {
            "folder": folder,
            "latlim": latlim,
            "lonlim": lonlim,
            "timelim": timelim,
            "product_name": product_name,
            "req_vars": req_vars,
        }

        dss[(source, product_name)] = dl_module.download(**args)

    if return_fps:
        for key, ds in dss.items():
            fp = ds.encoding["source"]
            ds = ds.close()
            dss[key] = fp

    return dss

def reverse_sources(sources):
    out = dict()
    for var, sublist in sources.items():
        if not isinstance(sublist, list):
            continue
        for source, product in sublist:
            if (source, product) in out.keys():
                out[(source, product)].append(var)
            else:
                out[(source, product)] = [var]
    return out

if __name__ == "__main__":

    import datetime

    sources = {

        "ndvi":         [("MODIS", "MOD13Q1.061"), ("MODIS", "MYD13Q1.061")],
        "r0":           [("MODIS", "MCD43A3.061")],
        # "lst":          [("MODIS", "MOD11A1.061"), ("MODIS", "MYD11A1.061")],

        "z":            [("SRTM", "30M")],
        "p":            [("CHIRPS", "P05")],
        "ra":           [("MERRA2", "M2T1NXRAD.5.12.4")],
        "t_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
        "t_air_max":    [("MERRA2", "M2I1NXASM.5.12.4")],
        # "u2m":          [("MERRA2", "M2I1NXASM.5.12.4")],
        # "v2m":          [("MERRA2", "M2I1NXASM.5.12.4")],
        # "qv":           [("MERRA2", "M2I1NXASM.5.12.4")],
        # "wv":           [("MERRA2", "M2I1NXASM.5.12.4")],
        # "p_air":        [("MERRA2", "M2I1NXASM.5.12.4")],
        # "p_air_0":      [("MERRA2", "M2I1NXASM.5.12.4")],

        # "lw_offset":    [("STATICS", "lw_offset")],
        # "lw_slope":     [("STATICS", "lw_slope")],
        # "r0_bare":      [("STATICS", "r0_bare")],
        # "r0_full":      [("STATICS", "r0_full")],
        # "rn_offset":    [("STATICS", "rn_offset")],
        # "rn_slope":     [("STATICS", "rn_slope")],
        # "t_amp_year":   [("STATICS", "t_amp_year")],
        # "t_opt":        [("STATICS", "t_opt")],
        # "vpd_slope":    [("STATICS", "vpd_slope")],
        # "z_oro":        [("STATICS", "z_oro")],

        # "level_name":   "sideloading",
    }

    folder = r"/Users/hmcoerver/Downloads/pywapor_test"
    # latlim = [26.9, 33.7]
    # lonlim = [25.2, 37.2]
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 7, 1), datetime.date(2020, 8, 1)]

    dss0 = collect_sources(folder, sources, latlim, lonlim, timelim)
