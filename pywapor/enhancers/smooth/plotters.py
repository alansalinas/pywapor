import matplotlib.pyplot as plt
import os
import warnings
import numpy as np
from matplotlib import ticker

def count_map(ax, ds):
    ds["ndvi"].count(dim = "time").rename("count").plot(ax = ax)
    ax.set_title("")

def cve_map(ax, ds):
    ds["cves"].sel(lmbda = ds["lmbda_sel"]).rename("cve").plot(ax = ax)
    ax.set_title("")

def plot_point(ax, point_ds, ylim = [-0.2, 1], t_idx = None, title = True, xtick = True):
    cmap = plt.cm.get_cmap('tab10')
    handles = []
    for i, sensor_name in enumerate(np.unique(point_ds["sensor"].values)):
        obj = ax.scatter(
                point_ds.time[point_ds.sensor == sensor_name], 
                point_ds.ndvi[point_ds.sensor == sensor_name], 
                color = cmap(i),
                label = sensor_name,
                )
        handles.append(obj)
    obj = ax.plot(point_ds["time"], point_ds["ndvi_smoothed"], label = "smoothed", color = cmap(i + 1))
    if not isinstance(ylim, type(None)):
        ax.set_ylim(ylim)
    if not isinstance(t_idx, type(None)):
        ax.plot([point_ds["time"][t_idx].values] * 2, ylim, color = "black", linewidth = 5, alpha = 0.3)
    handles.append(obj[0])
    ax.grid()
    ax.set_ylabel("ndvi [-]")
    ax.set_facecolor("lightgray")
    if not xtick:
        ax.set_xticklabels([])
        ax.set_xlabel("")
    else:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha="right")
        ax.set_xlabel("time [date]")
    if title:
        title_parts = []
        if "pname" in point_ds.attrs.keys():
            name = getattr(point_ds, "pname", "")
            title_parts.append(r"$\bf{" + name + r"}$")
        if "lmbda_sel" in point_ds.data_vars:
            lmbda = point_ds['lmbda_sel']
            title_parts.append(f"Lambda: {lmbda:.0E}")
        if "cves" in point_ds.data_vars and "lmbda_sel" in point_ds.data_vars:
            cve = point_ds['cves'].sel(lmbda = point_ds['lmbda_sel']).values
            title_parts.append(f"CVE: {cve:.4f}")
        if "x" in point_ds.coords:
            lat = float(point_ds.y.values)
            title_parts.append(f"Lat.: {lat:.3f}")
        if "y" in point_ds.coords:
            lon = float(point_ds.x.values)
            title_parts.append(f"Lon.: {lon:.3f}")
        full_title = ", ".join(title_parts)
        ax.set_title(full_title)
    ax.legend(handles = handles, ncols = 5, loc = "lower center")
    return handles

def plot_map(ax, da, points = None, ylim = [-1.0, 1.], ytick = True, xtick = True):
    assert da.ndim == 2
    im = ax.pcolormesh(da.x, da.y, da, cmap = "RdBu_r", vmin = ylim[0], vmax = ylim[1])
    ax.grid()
    ax.set_facecolor("lightgray")
    ax.set_title(da.name)
    if not isinstance(points, type(None)):
        ax.scatter(points[0], points[1], marker = "o", s = 100, edgecolors="black", linewidths = 2, c = [(1,1,1,0)], zorder = 100)
        for x,y,name in zip(*points):
            ax.annotate(name, (x,y), (30,0), textcoords = 'offset pixels', fontsize = 15)
    if not ytick:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel("Lat. [DD]")
    if not xtick:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel("Lon. [DD]")
    # ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=None))
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.3f}"))
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    return im

def plot_overview(ds, points, t_idx, folder):
    panels = [["A", "B"]]
    for i in range(len(points[0])):
        panels.append([str(i), str(i)])

    fig, axes = plt.subplot_mosaic(panels, figsize = (10, (len(points[0]) + 1)*4), dpi = 300)

    im1 = plot_map(axes["A"], ds["ndvi"].isel(time=t_idx), points=points)
    im2 = plot_map(axes["B"], ds["ndvi_smoothed"].isel(time=t_idx), ytick = False, points = points)

    for i, (lon, lat, name) in enumerate(zip(*points)):
        xtick = False if i+1 < len(points[0]) else True
        point_ds = ds.sel({"x": lon, "y": lat}, method = "nearest")
        point_ds = point_ds.assign_attrs({"pname": name})
        _ = plot_point(axes[str(i)], point_ds, t_idx = t_idx, title = True, xtick = xtick)

    fig.colorbar(im1)
    fig.colorbar(im2)
    fig.subplots_adjust(left=None, bottom=None, right=None, top=0.95, wspace=None, hspace=0.3)
    fig.suptitle(np.datetime_as_string(ds.time[t_idx], unit='m'))
    fig.savefig(os.path.join(folder, f"{t_idx:>03}.png"))
    plt.close(fig)

# ffmpeg -f image2 -framerate 5 -i %03d.png -vcodec libx264 -crf 22 video.mp4

# from PIL import Image
# import glob

# def make_gif(frame_folder, out_fh):
#     frames = [Image.open(image) for image in np.sort(glob.glob(frame_folder))]
#     frame_one = frames[0]
#     frame_one.save(out_fh, format="mp4", append_images=frames,
#                save_all=True, duration=int(len(frames)/2), loop=0)
    
# make_gif(r"/Users/hmcoerver/Local/test/*.png", r"/Users/hmcoerver/Local/test.mp4")
