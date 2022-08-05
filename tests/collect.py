#%%
import pywapor
import inspect
import datetime
import pandas as pd
import types
import functools
import csv
import os
import numpy as np

csv_lines = list()

workdir = r"/Users/hmcoerver/On My Mac/create_table"

# Determine output paths.
md_out = os.path.join(os.path.split(pywapor.__path__[0])[0], "docs", "mds", "data.md")
csv_out = os.path.join(os.path.split(pywapor.__path__[0])[0], "tests", "products_table.csv")

#%%

timelim1 = [datetime.date(2019, 7, 1), datetime.date(2019, 7, 15)]
timelim2 = [datetime.date(2019, 7, 1), datetime.date(2019, 7, 5)]
timelim3 = [datetime.date(2019, 7, 1), datetime.date(2019, 10, 1)]
timelim4 = [datetime.date(2022, 7, 1), datetime.date(2022, 7, 15)]
timelim5 = [datetime.date(2022, 7, 1), datetime.date(2022, 7, 3)]

latlim = [28.9, 29.7]
lonlim = [30.2, 31.2]

sources = {
    'GEOS5':        [pywapor.collect.product.GEOS5,     timelim2],
    'STATICS':      [pywapor.collect.product.STATICS,   timelim2],
    'MODIS':        [pywapor.collect.product.MODIS,     timelim3],
    'MERRA2':       [pywapor.collect.product.MERRA2,    timelim2],
    'GLOBCOVER':    [pywapor.collect.product.GLOBCOVER, timelim2],
    'CHIRPS':       [pywapor.collect.product.CHIRPS,    timelim2],
    'SRTM':         [pywapor.collect.product.SRTM,      timelim2],
    'PROBAV':       [pywapor.collect.product.PROBAV,    timelim1],
    'ERA5':         [pywapor.collect.product.ERA5,      timelim1],
    'SENTINEL2':    [pywapor.collect.product.SENTINEL2, timelim4],    
    'SENTINEL3':    [pywapor.collect.product.SENTINEL3, timelim5],
    'VIIRSL1':      [pywapor.collect.product.VIIRSL1,   timelim5]
}

overwrite_temp_res = {
    'sis-agrometeorological-indicators': "D"
}

# Loop over different sources.
for source, (mod, timelim) in sources.items():

    # Parse the `default_vars`-function code content.
    lines = inspect.getsourcelines(mod.default_vars)

    # Read the lines from `default_vars` to substract the availabel products and variables.
    idx1 = [i for i, x in enumerate(lines[0]) if "variables = " in x][0]
    idx2 = [i - 1 for i, x in enumerate(lines[0]) if "req_dl_vars = " in x][0]
    variables = None
    exec("".join(lines[0][idx1:idx2]).replace("  ", ""))

    # Read the lines from `default_vars` to subtract the `req_dl_vars`
    idx1 = [i for i, x in enumerate(lines[0]) if "req_dl_vars = " in x][0]
    idx2 = [i - 1 for i, x in enumerate(lines[0]) if "out = " in x][0]
    req_dl_vars = None
    exec("".join(lines[0][idx1:idx2]).replace("  ", ""))

    # List the available products for this `source`.
    products = variables.keys()

    # Loop over the products.
    for product_name in products:

        # List the available variables for this product.
        req_vars = list(req_dl_vars[product_name].keys())
        
        # List the necessary post-processors.
        pp = mod.default_post_processors(product_name, req_vars)

        print(source, product_name, req_vars)

        args = {
            "folder": workdir,
            "latlim": latlim,
            "lonlim": lonlim,
            "timelim": timelim, 
            "product_name": product_name,
            "req_vars": req_vars,
        }

        # Run the download function.
        ds = mod.download(**args)

        # Determine the spatial resolution of the dataset.
        sr = ds.rio.resolution()
        if np.isclose(abs(sr[0]), abs(sr[1]), atol = 0.00001):
            spat_resolution = f"{abs(sr[0]):.4f}°"
        else:
            spat_resolution = f"{abs(sr[0]):.4f}° x {abs(sr[1]):.4f}°"

        # Determine the temporal resolution.
        if "time" in ds.coords and ds.time.size >= 3:
            temp_resoltuion = pd.infer_freq(ds.time.values)
        else:
            temp_resoltuion = "N/A"
        if isinstance(temp_resoltuion, type(None)):
            temp_resoltuion = "N/A"
        if temp_resoltuion == "N/A" and product_name in overwrite_temp_res.keys():
            temp_resoltuion = overwrite_temp_res[product_name]

        # Loop over the downloaded variables.
        for var in req_vars:

            # Determine from which original variables the `var` is calculated.
            input_vars = ", ".join(req_dl_vars[product_name][var])

            # List the enhancers applied to this variable.
            enhance_names = list()
            for x in pp[var]:
                if isinstance(x, types.FunctionType):
                    enhance_names.append(f"`{x.__name__}`")
                elif isinstance(source, functools.partial):
                    enhance_names.append(f"`{x.func.__name__}`")
            enhancers = ", ".join(enhance_names)

            # Store everything in a list to be saved in a csv later.
            csv_line = [source, product_name, spat_resolution, temp_resoltuion, var, enhancers, input_vars]
            csv_lines.append(csv_line)

        print([x for x in ds.coords])
        print([x for x in ds.data_vars])
        print(ds.rio.crs)
        print("")

# Define headers.
header = ["source", "product", "spatial res.", "temp. res.", "variable", "processors", "inputs"]

# Write to csv.
with open(csv_out, 'w') as csvfile:
    filewriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    filewriter.writerow(header)
    for row in csv_lines:
        filewriter.writerow(row)

#%%

# Process the csv to markdown for publication in documentation.
df = pd.read_csv(csv_out, delimiter=";")
df = df.drop(["processors", "inputs"], axis = 1)
df = df.fillna("N/A")
df = df.groupby(["source", "product", "spatial res.", "temp. res."], as_index=False).agg({'variable': ', '.join,})
row = df[df.source == "STATICS"].groupby(["source", "spatial res.", "temp. res."], as_index=False).agg({'variable': ', '.join, 'product': ', '.join})
df = df[df.source != "STATICS"]
df = pd.concat([df, row], ignore_index=True).sort_values("source")

df.to_markdown(md_out.replace(".md", ".rst"), tablefmt="grid", index="never")
# %%
