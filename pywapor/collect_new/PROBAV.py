import os 
from pywapor.collect import accounts
import datetime
import pandas as pd
import pywapor.collect.accounts as accounts
from pywapor.general.logger import log
import re
import requests
from pywapor.collect_new.requests import download_urls, crawl

def download(folder, latlim, lonlim, timelim, product_name,
                variables = None, post_processors = None):
    
    dates = pd.date_range(timelim[0], timelim[1], freq="D")

    coords = [str(lonlim[0]), 
                str(latlim[0]), 
                str(lonlim[1]), 
                str(latlim[1])]

    base_url = f"https://www.vito-eodata.be/PDF/datapool/Free_Data/PROBA-V_100m/{product_name}"
    coord_req = "?coord=" + ",".join(coords)
    url = os.path.join(base_url, coord_req)

    session = requests.sessions.Session()
    session.auth = accounts.get("VITO")

    urls = find_tiles(url, dates, session)

    session = requests.sessions.Session()
    session.auth = accounts.get("VITO")
    fps = download_urls(urls, folder, session, parallel = 6)

    return fps

def find_tiles(url, dates, session):

    log.info("--> Searching PROBAV tiles.")

    regex = "https:.*\/\d{4}\/\?coord="
    filter_regex = "\d{4}(?=\/\?coord=)"
    urls = {"_": url}
    label_filter = dates.strftime("%Y")
    years = crawl(urls, regex, filter_regex, session, label_filter = label_filter)

    regex = "https:.*\/\d{4}\/\d{2}\/\?coord="
    filter_regex = "\d{4}\/\d{2}(?=\/\?coord=)"
    label_filter = dates.strftime("%Y%m")
    months = crawl(years, regex, filter_regex, session, label_filter = label_filter)

    regex = "https:.*\/\d{4}\/\d{2}\/\d{2}\/\?coord="
    filter_regex = "\d{4}\/\d{2}\/\d{2}(?=\/\?coord=)"
    label_filter = dates.strftime("%Y%m%d")
    days = crawl(months, regex, filter_regex, session, label_filter = label_filter)

    regex = "https:.*\/\d{4}\/\d{2}\/\d{2}\/.*\/\?coord="
    filter_regex = "\d{4}\/\d{2}\/\d{2}(?=\/.*\/\?coord=)"
    label_filter = dates.strftime("%Y%m%d")
    prods = crawl(days, regex, filter_regex, session, label_filter = label_filter)

    regex = ".*\.HDF5"
    filter_regex = "\d{8}(?=.*\.HDF5)"
    fns = crawl(prods, regex, filter_regex, session, list_out = True)

    dl_urls = [os.path.join(re.sub("\?coord=.*","",prods[date_str]), fn) for date_str, fn, in fns]

    log.info(f"--> Found {len(dl_urls)} PROBAV tiles.")

    return dl_urls

if __name__ == "__main__":

    product_name = "S5_TOC_100_m_C1"

    folder = r"/Users/hmcoerver/Downloads/merra2"

    latlim = [26.9, 33.7]
    lonlim = [25.2, 37.2]
    timelim = [datetime.date(2020, 7, 1), datetime.date(2020, 8, 15)]

    variables = None
    post_processors = None

    fps = download(folder, latlim, lonlim, timelim, product_name,
                variables = variables, post_processors = post_processors)


