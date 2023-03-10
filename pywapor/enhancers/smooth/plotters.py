import matplotlib.pyplot as plt
import os
import string
import warnings
import numpy as np
from matplotlib import ticker, colors
import glob

def plot_point(ax, point_ds, ylim = [-0.2, 1], t_idx = None, title = True, xtick = True):
    cmap = plt.cm.get_cmap('tab10')
    handles = []
    for i, sensor_name in enumerate(np.unique(point_ds["sensor"].values)):
        X = point_ds.time[point_ds.sensor == sensor_name]
        Y = point_ds.ndvi[point_ds.sensor == sensor_name]
        if Y.count("time").values > 0:
            obj = ax.scatter(
                    X, Y,
                    color = cmap(i),
                    label = point_ds.sensor.attrs[str(sensor_name)],
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
            title_parts.append(f"Lambda: {lmbda.values:.0E}")
        if "cves" in point_ds.data_vars and "lmbda_sel" in point_ds.data_vars:
            cve = point_ds['cves'].sel(lmbda = point_ds['lmbda_sel'], method = "nearest").values
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

def plot_map(ax, da, points = None, cmap = "RdBu_r", ylim = [-1.0, 1.0], ytick = True, xtick = True, norm = None):
    assert da.ndim == 2
    im = ax.pcolormesh(da.x, da.y, da, cmap = cmap, vmin = ylim[0], vmax = ylim[1], norm=norm)
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
        ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.3f}"))
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    if not xtick:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel("Lon. [DD]")
        ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.3f}"))
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))

    return im

def plot_overview(ds, points, t_idx, folder):

    if not os.path.exists(folder):
        os.makedirs(folder)

    panels = [["A", "B"]]
    for i in range(len(points[0])):
        panels.append([str(i), str(i)])

    fig, axes = plt.subplot_mosaic(panels, figsize = (10, (len(points[0]) + 1)*4), dpi = 16**2)

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
    fig.savefig(os.path.join(folder, f"{t_idx:>06}.png"))
    plt.close(fig)


def plot_stats(ds, folder):

    if not os.path.exists(folder):
        os.makedirs(folder)

    npanels = np.sum(np.isin(["lmbda_sel", "ndvi", "cves"], ds.data_vars))
    panels = list(string.ascii_uppercase[:npanels])

    fig, axes = plt.subplot_mosaic([panels], figsize = (npanels * 5, 4), dpi = 300)
    
    if "ndvi" in ds.data_vars:
        ds["counts"] = ds["ndvi"].count(dim = "time")
        im1 = plot_map(axes[panels.pop(0)], ds["counts"], ylim = [None, None], cmap = "viridis")
        fig.colorbar(im1, label = "No. measurements [-]")

    if np.all(np.isin(["cves", "lmbda_sel"], ds.data_vars)):
        if "lmbda" in ds["cves"].dims:
            ds["cve"] = ds["cves"].sel(lmbda = ds["lmbda_sel"], method = "nearest")
        else:
            ds["cve"] = ds["cves"]
        im2 = plot_map(axes["B"], ds["cve"], ylim = [None, None], cmap = "viridis", ytick = False)
        fig.colorbar(im2, label = "Cross-val. Standard Error [-]")
    
    if "lmbda_sel" in ds.data_vars:
        im3 = plot_map(axes["C"], ds["lmbda_sel"], ylim = [None, None], cmap = "viridis", ytick = False, norm = colors.LogNorm())
        fig.colorbar(im3, label = "Lambda [-]")

    fig.savefig(os.path.join(folder, f"stats.png"))
    plt.close(fig)

def create_video(files, video_fh, fps = 5, remove_files = True):

    try:
        import imageio.v2 as imageio
    except ImportError:
        print("--> Install `imageio` to automatically create a video")
        return

    try:
        with imageio.get_writer(video_fh, fps = fps) as writer:
            for im in files:
                writer.append_data(imageio.imread(im))
        if remove_files:
            for fn in files:
                try:
                    os.remove(fn)
                except PermissionError:
                    continue
    except ValueError as e:
        msg = getattr(e, "args", ["None"])
        if "Based on the extension, the following" in msg[0]:
            fn, ext = os.path.splitext(video_fh)
            if ext == ".mp4":
                print("--> Creating `.gif` file, install `imageio_ffmpeg` to create `.mp4`.")
                create_video(files, fn + ".gif", fps = fps)
            else:
                print("--> Unable to create video with requested extension.")
                print(e)
        else:
            raise e

# ffmpeg -f image2 -framerate 5 -i %03d.png -vcodec libx264 -crf 22 video.mp4

# from PIL import Image
# import glob

# def make_gif(frame_folder, out_fh):
#     frames = [Image.open(image) for image in np.sort(glob.glob(frame_folder))]
#     frame_one = frames[0]
#     frame_one.save(out_fh, format="mp4", append_images=frames,
#                save_all=True, duration=int(len(frames)/2), loop=0)
    
# make_gif(r"/Users/hmcoerver/Local/test/*.png", r"/Users/hmcoerver/Local/test.mp4")

# fileList = []
# for file in os.listdir(path):
#     if file.startswith(name):
#         complete_path = path + file
#         fileList.append(complete_path)

# writer = imageio.get_writer('test.mp4', fps=20)

# for im in fileList:
#     writer.append_data(imageio.imread(im))
# writer.close()