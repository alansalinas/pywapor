import pywapor
import os
import numpy as np
import glob
import warnings
import json
from functools import partial
from pywapor.general.logger import log, adjust_logger
from pywapor.collect.downloader import collect_sources
from pywapor.general.levels import pre_et_look_levels
from pywapor.general.compositer import time_bins
from pywapor.general.processing_functions import adjust_timelim_dtype, func_from_string, open_ds

class Configuration():

    variable_categories = {
        "optical": ["ndvi", "r0", "mndwi", "psri", "vari_red_edge", "bsi", "nmdi", "green", "nir"],
        "thermal": ["bt", "lst"],
        "solar radiation": ["ra"],
        "precipitation": ["p"],
        "elevation": ["z", "slope", "aspect"],
        "meteorological": ["t_air", "t_air_min", "t_air_max", "u", "vp", "p_air", "p_air_0", "u2m", "v2m", "qv", "wv", "t_dew"],
        "statics": ["lw_slope", "lw_offset", "z_obst_max", "rs_min", "land_mask", "vpd_slope", "t_opt", "t_amp", "t_amp_year", "rn_slope", "rn_offset", "z_oro"],
        "soil moisture": ["se_root"],
    }

    category_variables = {}
    for x in [{var: cat for var in varis} for cat, varis in variable_categories.items()]:
        category_variables = {**category_variables, **x}

    se_root_vars = [
        [("ndvi",)], 
        [("bt", "lst"), ("lst",), ("bt",)], 
        [("u",), ("u2m", "v2m")],
        [("p_air",)], 
        [("p_air_0",)], 
        [("t_air",)], 
        [("wv",)], 
        [("t_dew",), ("qv",)]
        ]
    
    et_look_vars = [
        [("ndvi",)],
        [("r0",)],

        [("se_root",)], 

        [("p",)],
        [("z",)],
        [("ra",)],

        [("u",), ("u2m", "v2m")],
        [("p_air",)], 
        [("p_air_0",)], 
        [("t_air",)],
        [("t_air_min",)], 
        [("t_air_max",)], 
        [("t_dew",), ("qv",)],

        [("rn_offset",)],
        [("rn_slope",)],
        [("t_amp",),],
        [("t_opt",)],
        [("vpd_slope",)],
        [("lw_offset",)],
        [("lw_slope",)],

        [("rs_min",), ()],
        [("land_mask",), ()],
        [("z_obst_max",), ()],
        [("z_oro",), ()],
        ]
    
    @staticmethod
    def source_func(x):
        return x if "FILE:" in x else x.split(".")[0]
    
    @staticmethod
    def pname_func(x):
        return "from_file" if "FILE:" in x else ".".join(x.split(".")[1:])

    def __init__(self, full = None, summary = None, se_root = None, 
                 et_look = None):
        self.full = full
        self.summary = summary
        self.se_root = se_root
        self.et_look = et_look

    def __repr__(self):
        summary = self.summary.copy()

        example_product = summary.pop('_EXAMPLE_')
        whittaker = summary.pop('_WHITTAKER_')
        sharpen = summary.pop('_ENHANCE_')

        base = "\n".join([f"--> {cat.title()} data:\n    > `{'`, `'.join(sources)}`" for cat, sources in summary.items()])

        additional = f"""
--> Example Product:
    > {example_product}
--> Applying whittaker to variables from sources:
    > `{'` and `'.join(whittaker)}`
--> Sharpening these variables:
    > `{'` and `'.join(sharpen)}`"""
        return base + additional

    @staticmethod
    def default_enhancers(prod, var):
        source = Configuration.source_func(prod)
        if "FILE:" in source:
            return []
        product_name = Configuration.pname_func(prod)
        mod = getattr(pywapor.collect.product, source)
        x = mod.default_post_processors(product_name, [var])[var]

        funcs = list()
        for x_ in x:
            if isinstance(x_, partial):
                funcs.append({
                    "func": f"{x_.func.__module__}.{x_.func.__name__}",
                    "args": x_.args,
                    "keywords": x_.keywords,
                })
            else:
                funcs.append({"func": f"{x_.__module__}.{x_.__name__}"})
        
        return funcs

    @classmethod
    def from_name(cls, name):
        log.info(f"--> Searching configuration for `{name}`.").add()

        folder = os.path.realpath(pywapor.__path__[0])
        available_configs = glob.glob(os.path.join(folder, "configs", "*.json"))
        available_names = sorted([os.path.splitext(os.path.split(x)[-1])[0] for x in available_configs])
        log.info(f"> Available configurations are `{'`, `'.join(available_names)}`")
    
        synonyms = {
            "level_1": "WaPOR2_level_1",
            "level_2": "WaPOR2_level_2",
            "level_3": "WaPOR2_level_3",
            "level_1_v2": "WaPOR2_level_1",
            "level_2_v2": "WaPOR2_level_2",
            "level_3_v2": "WaPOR2_level_3",
            "level_1_v3": "WaPOR3_level_1",
            "level_2_v3": "WaPOR3_level_2",
            "level_3_v3": "WaPOR3_level_3",
        }

        name_ = synonyms.get(name, name)

        selection = None
        for fh in available_configs:
            fn = os.path.splitext(os.path.split(fh)[-1])[0]
            if fn == name_:
                selection = fh
                break

        if isinstance(selection, type(None)):
            raise ValueError(f"No configuration found for name `{name}`, choose one from `{'`, `'.join(available_names)}`.")
        else:
            config = Configuration.from_json(selection)

        log.sub().info("--> Configuration set.")

        return config
    
    def to_json(self, fh):
        out = {
            "summary": self.summary, 
            "full": self.full, 
            "se_root": self.se_root, 
            "et_look": self.et_look,
        }
        class EncodeSetTuple(json.JSONEncoder):
            def encode(self, obj):
                def hint_tuples(item):
                    if isinstance(item, tuple):
                        return {'type': "tuple", 'value': list(item)}
                    if isinstance(item, list):
                        return [hint_tuples(e) for e in item]
                    if isinstance(item, set):
                        return {'type': "set", 'value': list(item)}
                    if isinstance(item, dict):
                        return {key: hint_tuples(value) for key, value in item.items()}
                    else:
                        return item
                return super(EncodeSetTuple, self).encode(hint_tuples(obj))
        
        json_string = json.dumps(out, cls = EncodeSetTuple)
        with open(fh, "w") as x:
            x.write(json_string)

    @classmethod
    def from_json(cls, fh):
        log.info(f"--> Loading configuration from `{os.path.split(fh)[1]}`.").add()
        def decode_set(dct):
            if dct.get("type", None) == "set":
                return set(dct.get('value'))
            if dct.get("type", None) == "tuple":
                return tuple(dct.get('value'))
            return dct
        with open(fh, "r") as fp:
            config_ = json.load(fp, object_hook=decode_set)
        
        config = cls(full = config_["full"], summary = config_["summary"], se_root = config_["se_root"], 
                 et_look = config_["et_look"])

        config.validate()
        log.sub().info("--> Configuration loaded.")
        return config
    
    @classmethod
    def from_summary(cls, summary):
        log.info(f"--> Creating configuration from summary.").add()
        example_product = summary.pop('_EXAMPLE_', '')
        temporal_interp = summary.pop('_TEMP_INTERP_', {})
        enhance = summary.pop('_ENHANCE_', {})

        full = dict()

        unique_enhancers = set([item for row in list(enhance.values()) for item in row])
        extra_vars = list()
        for enhancer in unique_enhancers:
            func = func_from_string(enhancer)
            extra_vars += [[(x,)] for x in func.__kwdefaults__["req_vars"]]

        for var_group in cls.et_look_vars + cls.se_root_vars + extra_vars:
            for var_ in var_group:
                for var in var_:
                    cat = cls.category_variables[var]
                    possible_prods = summary[cat]
                    products = []
                    t_interps = list()
                    for prod in possible_prods:
                        valid = cls.has_var(prod, var)
                        if valid:

                            enhancers = cls.default_enhancers(prod, var)

                            products += [{
                                "source": cls.source_func(prod), 
                                "product_name": cls.pname_func(prod), 
                                "enhancers": enhancers,
                                "is_example": prod == example_product,
                                }]
                        t_interps += [temporal_interp.get(prod, [])]
                    
                    if len(products) == 0:
                        continue

                    if len(t_interps) == 0:
                        temporal_interp_ = "linear"
                    elif len(t_interps) == 1:
                        temporal_interp_ = t_interps[0]
                    else:
                        log.warning("")
                        temporal_interp_ = t_interps[0]

                    variable_enhancers = enhance.get(var, [])

                    composite_type = "mean" if var.split("_")[-1] not in ["min", "max"] else var.split("_")[-1]

                    full[var] = {
                        "products": products,
                        "temporal_interp": temporal_interp_,
                        "variable_enhancers": variable_enhancers,
                        "spatial_interp": "bilinear",
                        "composite_type": composite_type,
                    }

        config = cls(full = full)
        config.validate()
        config.summarize()
        config.update_se_root_config()
        config.update_et_look_config()
        log.sub().info("--> Configuration created.")
        return config

    @staticmethod
    def has_var(prod, var, verbose = True):
        source = Configuration.source_func(prod)
        product_name = Configuration.pname_func(prod)
        if "FILE:" in source:
            if not verbose:
                log.info(f"> Variable `{var}` will be loaded from a file `{source}`.")
            return True
        mod = getattr(pywapor.collect.product, source)
        try:
            mod.default_vars(product_name, [var])
            valid = True
        except TypeError:
            valid = False
            if not verbose:
                log.warning(f"> {source}.{product_name} does not have a variable called `{var}`.") 
        return valid
    
    def validate(self):
        log.info("--> Validating configuration.").add()
        valids = []
        for var in self.full:
            for product in self.full[var]["products"]:
                prod = f"{product['source']}.{product['product_name']}"
                valid = Configuration.has_var(prod, var, verbose = False)
                valids.append(valid)
        if all(valids):
            log.info("> All specified variables sources are valid.")
        log.sub()
        return valid

    def update_se_root_config(self):
        log.info("--> Making configuration for SE_ROOT.").add()
        
        sharpened_vars = {var for var, config in self.full.items() if any(["sharpen" in x for x in config["variable_enhancers"]])}

        sharpeners = ["mndwi", "psri", "vari_red_edge", "bsi", "nmdi", "green", "nir"]
    
        se_root_config = dict()
        for var_group in self.se_root_vars:
            for var in var_group:
                if all([var_ in self.full.keys() for var_ in var]):
                    for var_ in var:
                        se_root_config[var_] = self.full[var_].copy()
                        _ = se_root_config[var_].pop("composite_type")
                    break
                elif var == var_group[-1]:
                    substring = [f"[`{'` and `'.join(x)}`]" for x in var_group]
                    log.warning(f"--> No configuration found for {' or '.join(substring)}.")
                else:
                    continue
        
        to_sharpen = set(sharpened_vars) & set(se_root_config)
        if to_sharpen:
            for var in sharpeners:
                if var in self.full.keys():
                    se_root_config[var] = self.full[var]
                else:
                    log.warning(f"--> No configuration found for `{var}` (required for sharpening of `{'` and `'.join(to_sharpen)}`).")

        self.se_root = se_root_config
        log.sub()

    def update_et_look_config(self):
        log.info("--> Making configuration for ET_LOOK.").add()
        
        sharpened_vars = {var for var, config in self.full.items() if any(["sharpen" in x for x in config["variable_enhancers"]])}

        sharpeners = ["mndwi", "psri", "vari_red_edge", "bsi", "nmdi", "green", "nir"]
    
        et_look_config = dict()
        for var_group in self.et_look_vars:
            for var in var_group:
                if all([var_ in self.full.keys() for var_ in var]):
                    for var_ in var:
                        et_look_config[var_] = self.full[var_].copy()
                    if len(var) == 0 and var == var_group[-1]:
                        substring = [f"[`{'` and `'.join(x)}`]" for x in var_group[:-1]]
                        log.info(f"--> No configuration found for {' or '.join(substring)}, will use constant value.")
                    break
                elif var == var_group[-1]:
                    substring = [f"[`{'` and `'.join(x)}`]" for x in var_group]
                    log.warning(f"--> No configuration found for {' or '.join(substring)}.")
                else:
                    continue
        
        to_sharpen = set(sharpened_vars) & set(et_look_config)
        if to_sharpen:
            for var in sharpeners:
                if var in self.full.keys():
                    et_look_config[var] = self.full[var]
                else:
                    log.warning(f"--> No configuration found for `{var}` (required for sharpening of `{'` and `'.join(to_sharpen)}`).")

        self.et_look = et_look_config
        log.sub()

    def summarize(self):
        log.info("--> Making summary of configuration.")
        
        summary = dict()
        for var, config in self.full.items():
            sources = set([f'{x["source"]}.{x["product_name"]}' for x in config["products"]])
            cat = self.category_variables.get(var, "unknown")
            if cat == "other":
                log.warning(f"--> Unknown variable `{var}` found.")
            if cat in list(summary.keys()):
                summary[cat].union(sources)
            else:
                summary[cat] = sources

        whittaker = dict()
        sharpen = set()
        for k,v in self.full.items():
            example_ = [f"{prod['source']}.{prod['product_name']}" for prod in v["products"] if prod["is_example"]]
            if len(example_) > 0:
                example = example_[0]
            if isinstance(v["temporal_interp"], dict):
                if v["temporal_interp"].get("method", "") == "whittaker":
                    for prod in v["products"]:
                        whittaker[f"{prod['source']}.{prod['product_name']}"] = v["temporal_interp"]
            if any(["sharpen" in x for x in v["variable_enhancers"]]):
                sharpen.add(k)

        summary["_WHITTAKER_"] = whittaker
        summary["_ENHANCE_"] = sharpen
        summary["_EXAMPLE_"] = example

        self.summary = summary

class Project():

    def __init__(self, project_folder, bb, period, configuration = None):

        assert bb[::2][0] < bb[::2][1], "Invalid Bounding-Box"
        self.lonlim = bb[::2]
        assert bb[1::2][0] < bb[1::2][1], "Invalid Bounding-Box"
        self.latlim = bb[1::2]
        if not os.path.isdir(project_folder):
            os.makedirs(project_folder)
        self.folder = project_folder
        self.period = adjust_timelim_dtype(period)
        self.configuration = configuration
        self.dss = None

        warnings.filterwarnings("ignore", message="invalid value encountered in power")

        adjust_logger(True, self.folder, "INFO")

        if os.path.isfile(os.path.join(self.folder, "se_root_in.nc")):
            self.se_root_in = open_ds(os.path.join(self.folder, "se_root_in.nc"))
        else:
            self.se_root_in = None

        se_root_outs = glob.glob(os.path.join(self.folder, "se_root_out*.nc"))
        if se_root_outs:
            self.se_root_out = open_ds(max(se_root_outs, key=os.path.getmtime))
        else:
            self.se_root_out = None

        if os.path.isfile(os.path.join(self.folder, "et_look_in.nc")):
            self.et_look_in = open_ds(os.path.join(self.folder, "et_look_in.nc"))
        else:
            self.et_look_in = None

        et_look_outs = glob.glob(os.path.join(self.folder, "et_look_out*.nc"))
        if et_look_outs:
            self.et_look_out = open_ds(max(et_look_outs, key=os.path.getmtime))
        else:
            self.et_look_out = None

        log.info("> PROJECT").add()
        log.info(self.__repr__())
        log.sub().info("< PROJECT")

    def __repr__(self):
        project_str = f"""--> Project Folder:
        > {self.folder}
    --> Period:
        > {self.period[0]} - {self.period[1]}
    --> Bounding-Box:\n
                 {self.latlim[1]:8.4f}
                 ┌─────────┐
                 │         │
       {self.lonlim[0]:9.4f} │         │{self.lonlim[1]:9.4f}
                 │         │
                 └─────────┘
                 {self.latlim[0]:8.4f}\n
    --> Configuration:
        > {self.configuration}
    --> pyWaPOR Version:
        > {pywapor.__version__}"""
        return project_str

    def load_configuration(self, name = None, summary = None, json = None):
        if not isinstance(name, type(None)) and not isinstance(summary, type(None)):
            raise ValueError("Only one of `level` and `summary` can be specified.")
        log.info("> CONFIGURATION").add()
        if not isinstance(name, type(None)):
            self.configuration = Configuration.from_name(name)
        elif not isinstance(summary, type(None)):
            summary = summary.copy()
            self.configuration = Configuration.from_summary(summary)
        elif not isinstance(json, type(None)):
            self.configuration = Configuration.from_json(json)
        else:
            raise ValueError("At least one of `level` or `summary` needs to be specified.")
        
        log.info("--> Summary:").add()
        summary_str = self.configuration.__repr__()

        for line in summary_str.split("\n"):
            log.info(line)

        log.sub().sub().info("< CONFIGURATION")
        return self.configuration
    
    def set_passwords(self):

        log.info("> PASSWORDS").add()

        req_accounts = {
            "MODIS": "NASA",
            "SRTM": "NASA",
            "CHIRPS": "NASA",
            "MERRA2": "NASA",
            "TERRA": "TERRA",
            "ERA5": "ECMWF",
            "LANDSAT": "EARTHEXPLORER",
            "SENTINEL2": "COPERNICUS_DATA_SPACE",
            "SENTINEL3": "COPERNICUS_DATA_SPACE",
        }

        all_accounts = list()
        for v in self.configuration.summary.values():
            if isinstance(v, set):
                for source_product in v:
                    name = Configuration.source_func(source_product)
                    account = req_accounts.get(name, None)
                    if not isinstance(account, type(None)):
                        all_accounts.append(account)

        x = set(all_accounts)
        log.info(f"--> Accounts needed for `{'`, `'.join(x)}`.").add()

        for account in x:
            _ = pywapor.collect.accounts.get(account)

        log.sub().info("--> All set!")
        log.sub().info("< PASSWORDS")

    def download_data(self, buffer_timelim = True):
        assert not isinstance(self.configuration, type(None)), "Please load a configuration before continueing."

        log.info("> DOWNLOADER").add()

        if buffer_timelim:
            bins = time_bins(self.period, 1)
            adjusted_timelim = [bins[0], bins[-1]]
            buffered_timelim = [adjusted_timelim[0] - np.timedelta64(3, "D"), 
                                adjusted_timelim[1] + np.timedelta64(3, "D")]
        else:
            adjusted_timelim = self.period
            buffered_timelim = self.period

        self.dss, _ = collect_sources(
                            self.folder,
                            self.configuration.full,
                            self.latlim,
                            self.lonlim,
                            buffered_timelim,
                            )
        
        log.sub().info("< DOWNLOADER")
    
        return self.dss
    
    def run_pre_se_root(self, forced = False):
        if isinstance(self.se_root_in, type(None)) or forced:
            self.se_root_in = pywapor.pre_se_root.main(self.folder, self.latlim, self.lonlim, self.period, sources = self.configuration.se_root)
        else:
            log.info("> PRE_SE_ROOT").add()
            log.info(f"--> Re-using `{os.path.split(self.se_root_in.encoding['source'])[-1]}`.")
            log.sub().info("< PRE_SE_ROOT")
        return self.se_root_in

    def run_se_root(self, se_root_version = "v3"):
        self.se_root_out = pywapor.se_root.main(self.se_root_in, se_root_version = se_root_version)
        return self.se_root_out
    
    def run_pre_et_look(self, forced = False):
        if isinstance(self.et_look_in, type(None)) or forced:
            self.et_look_in = pywapor.pre_et_look.main(self.folder, self.latlim, self.lonlim, self.period, sources = self.configuration.et_look)
        else:
            log.info("> PRE_ET_LOOK").add()
            log.info(f"--> Re-using `{os.path.split(self.et_look_in.encoding['source'])[-1]}`.")
            log.sub().info("< PRE_ET_LOOK")
        return self.et_look_in

    def run_et_look(self, et_look_version = "v3"):
        self.et_look_out = pywapor.et_look.main(self.et_look_in, et_look_version = et_look_version)
        return self.et_look_out

if __name__ == "__main__":

    timelim = ["2022-01-01", "2023-12-31"]
    latlim = [21.9692194682626933, 21.9939120838340507]
    lonlim = [91.9371349243682801, 91.9657566608824339]
    bb = [91.9371349243682801, 21.9692194682626933, 91.9657566608824339, 21.9939120838340507]
    project_folder = r"/Users/hmcoerver/Local/pywapor_bgd"

    adjust_logger(True, project_folder, "INFO")

    all_configs = {
        "WaPOR2_level_1",
        "WaPOR2_level_2",
        "WaPOR2_level_3",
        "WaPOR2_level_1",
        "WaPOR2_level_2",
        "WaPOR2_level_3",
        "WaPOR3_level_2",
        "WaPOR3_level_3",
        "nrt",
        "all_in",
    }

    project = pywapor.Project(project_folder, bb, timelim)

    # project.load_configuration(name = "WaPOR3_level_2")

    summary = {
        # Define which products to use.
        'elevation': {'COPERNICUS.GLO30'},
        'meteorological': {'GEOS5.tavg1_2d_slv_Nx'},
        'optical': {'SENTINEL2.S2MSI2A_R20m'},
        'precipitation': {'CHIRPS.P05'},
        'solar radiation': {'ERA5.sis-agrometeorological-indicators'},
        'statics': {'STATICS.WaPOR3'},
        'thermal': {'VIIRSL1.VNP02IMG'},
        'soil moisture': {'FILE:{folder}{sep}se_root_out*.nc.from_file'},

        # Define which product to reproject the other products to.
        '_EXAMPLE_': 'SENTINEL2.S2MSI2A_R20m', 

        # Define any special functions to apply to a specific variable.
        '_ENHANCE_': {"bt": ["pywapor.enhancers.dms.thermal_sharpener.sharpen"],},

        # Choose which products should be gapfilled.
        '_WHITTAKER_': {
            'SENTINEL2.S2MSI2A_R20m': {'lmbdas': 1000.0, 'method': 'whittaker'}, 
            'VIIRSL1.VNP02IMG': {'a': 0.85, 'lmbdas': 1000.0, 'method': 'whittaker'}},
        }
    
    project.load_configuration(summary = summary)



    project.set_passwords()
    dss = project.download_data()
    
    se_root_in = project.run_pre_se_root()
    # se_root = project.run_se_root()
    # et_look_in = project.run_pre_et_look()
    # et_look = project.run_et_look()





