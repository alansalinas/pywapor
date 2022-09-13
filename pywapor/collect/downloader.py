from pywapor.general.logger import log, adjust_logger
import types
import functools
import numpy as np
import sys

def collect_sources(folder, sources, latlim, lonlim, timelim, return_fps = True):
    
    reversed_sources, reversed_enhancers = reverse_sources(sources)

    dss = dict()

    max_attempts = 3
    attempts = {k: 0 for k in reversed_sources.keys()}

    while not np.all([v >= max_attempts for v in attempts.values()]):

        reversed_sources = {k: v for k, v in reversed_sources.items() if attempts[k] < max_attempts}

        for (source, product_name), req_vars in reversed_sources.items():
            
            if isinstance(source, str):
                dl_module = __import__(f"pywapor.collect.product.{source}", 
                                    fromlist=[source])
                dler = dl_module.download
                log.info(f"--> Collecting `{'`, `'.join(req_vars)}` from {source}.{product_name}.")
                source_name = source
            elif isinstance(source, types.FunctionType):
                dler = source
                log.info(f"--> Collecting `{'`, `'.join(req_vars)}` from `{dler.__name__}`.")
                source_name = source.__name__
            elif isinstance(source, functools.partial):
                dler = source
                log.info(f"--> Collecting `{'`, `'.join(req_vars)}` from `{dler.func.__name__}`.")
                source_name = source.func.__name__

            log.add()
            
            args = {
                "folder": folder,
                "latlim": latlim,
                "lonlim": lonlim,
                "timelim": timelim,
                "product_name": product_name,
                "req_vars": req_vars,
                "post_processors": reversed_enhancers[(source, product_name)]
            }

            try:
                x = dler(**args)
                attempts[(source, product_name)] += max_attempts * 10
            except Exception as e:
                if attempts[(source, product_name)] < max_attempts - 1:
                    log.info(f"--> Attempt {attempts[(source, product_name)] + 1} of {max_attempts} failed, trying again after other sources have been collected. ({type(e).__name__}: {e}).")
                else:
                    log.info(f"--> Attempt {attempts[(source, product_name)] + 1} of {max_attempts} failed, giving up now. ({type(e).__name__}: {e}).")
                attempts[(source, product_name)] += 1
            else:
                if "time" in x.coords:
                    stime = np.datetime_as_string(x.time.values[0], unit = "m")
                    etime = np.datetime_as_string(x.time.values[-1], unit = "m")
                    log.add().info(f"> timesize: {x.time.size} [{stime}, ..., {etime}]").sub()
                dss[(source_name, product_name)] = x
            finally:
                log.sub()
                continue

    if return_fps:
        for key, ds in dss.items():
            fp = ds.encoding["source"]
            ds = ds.close()
            dss[key] = fp

    return dss

def reverse_sources(sources):
    reversed_sources = dict()
    reversed_enhancers = dict()
    for var, value in sources.items():
        for src in value["products"]:
            key = (src["source"], src["product_name"])
            enhancers = src["enhancers"]

            if key in reversed_sources.keys():
                reversed_sources[key].append(var)
                reversed_enhancers[key][var] = enhancers
            else:
                reversed_sources[key] = [var]
                reversed_enhancers[key] = {var: enhancers}

    return reversed_sources, reversed_enhancers

if __name__ == "__main__":

    import datetime

    sources = {

        "r0":           {"products": [{"source": "MODIS", "product_name": "MCD43A3.061", "enhancers": "default"}]},
        "z":            {"products": [{"source": "SRTM", "product_name": "30M", "enhancers": "default"}]},
        "p":            {"products": [{"source": "CHIRPS", "product_name": "P05", "enhancers": "default"}]},
        "lw_offset":    {"products": [{"source": "STATICS", "product_name": "WaPOR2", "enhancers": "default"}]},
    }

    folder = r"/Users/hmcoerver/Local/supersafe_dler"
    latlim = [28.9, 29.7]
    lonlim = [30.2, 31.2]
    timelim = [datetime.date(2020, 7, 1), datetime.date(2020, 7, 15)]

    adjust_logger(True, folder, "INFO")

    dss0 = collect_sources(folder, sources, latlim, lonlim, timelim)
