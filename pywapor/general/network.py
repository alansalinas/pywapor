from pyvis.network import Network
import json

def make_node(ds, var):

    if var in list(ds.variables):
        if "units" in ds[var].attrs.keys():
            units = ds[var].units
        else:
            units = ""

        if "long_name" in ds[var].attrs.keys():
            long_name = ds[var].long_name
        else:
            long_name = ""

        if "et_look_module" in ds[var].attrs.keys():
            group = ds[var].et_look_module
            shape = "dot"
            module = ds[var].et_look_module
        else:
            if len(ds[var].coords) == 0:
                group = "scalar"
            elif len(ds[var].coords) == 2:
                group = "temporal-constant"
            elif len(ds[var].coords) == 3:
                group = "temporal-input"
            elif len(ds[var].coords) > 3:
                group = "main-input"
            else:
                group = "unknown"
            shape = "square"
            module = "pre_et_look"

        if "calculated_with" in ds[var].attrs.keys():

            upstreams = ds[var].calculated_with
            if not isinstance(upstreams, list):
                upstreams = [upstreams]

            f = " = f(" + ", ".join(upstreams) + ")"
        else:
            f = ""


    else:
        long_name = ""
        units = ""
        module = ""
        group = ""
        module = ""
        f = ""
        shape = "dot"
        
    title = f"""
    <u>{long_name}</u><br><br>
    <b>{var}</b>{f} <br><br>
    Units: [{units}] <br>
    Module: {module}
    """

    return {"n_id": var, "label": var, "title": title, "group": group, "shape": shape}

def create_network(ds, fh, exclude = ["scalar", "temporal-input"]):
    """_summary_

    [
            "scalar", 
            # "temporal-constant", 
            "temporal-input",
            ]
    Parameters
    ----------
    fh : _type_
        _description_
    exclude : list, optional
        _description_, by default ["scalar", # "temporal-constant", "temporal-input"]
    """
    net = Network(width = "100%", height = "100%", layout = False, 
                    bgcolor='#E9ECF5', font_color='black', directed = True)

    options = {
        "physics": {
            "hierarchicalRepulsion": {
                "centralGravity": 0,
            },
        "minVelocity": 0.75,
        "solver": "hierarchicalRepulsion",
        },
        "edges": {
            "color": {
                "color": "#B5BFDE",
                "highlight": "black",
            },
            "selectionWidth": 4,
            "smooth": {
                "type": "cubicBezier",
            }
        },
        "layout": {
            "randomSeed": 460885,
            "improvedLayout": True,
        },
        "groups": {
            "main-input": {
                "color": {
                    "background": "#122332", 
                    "border": "#14293B", 
                    "highlight": {
                        "border": "#2A4257", 
                        "background": "#183751"
                        }, 
                    "hover": {
                        "background": "#2A4257", 
                        "border": "#183751"}
                        }, 
                "borderWidth": 1, 
                "borderWidthSelected": 2
                },
        },
        "interaction": {
            "hover": True,
            "tooltipDelay": 10,
            # "zoomSpeed": 3, doesnt work.
            "navigationButtons": False,
        }
    }

    net.set_options(json.dumps(options))

    for var in list(ds.variables):

        if "calculated_with" in ds[var].attrs.keys():

            upstreams = ds[var].calculated_with
            if not isinstance(upstreams, list):
                upstreams = [upstreams]

            if not var in net.get_nodes():
                kwargs = make_node(ds, var)
                if kwargs["group"] in exclude:
                    continue
                net.add_node(**kwargs)

            for upstream in upstreams:
                if not upstream in net.get_nodes():
                    kwargs = make_node(ds, upstream)
                    if kwargs["group"] in exclude:
                        continue
                    net.add_node(**kwargs)

                net.add_edge(upstream, var)

    net.show(fh)