import os
import urllib
from bs4 import BeautifulSoup
import urllib
from joblib import Parallel, delayed
from functools import partial
import re
import numpy as np
import tqdm
from requests.exceptions import HTTPError
import time
import socket
import glob
import requests
from pywapor.general.logger import log, adjust_logger
from cachetools import cached, TTLCache

@cached(cache=TTLCache(maxsize=2048, ttl=3600))
def find_paths(url, regex, node_type = "a", tag = "href", filter = None, session = None):
    if isinstance(session, type(None)):
        f = urllib.request.urlopen(url)
    else:
        file_object = session.get(url, stream = True)
        file_object.raise_for_status()
        f = file_object.content
    soup = BeautifulSoup(f, "lxml")
    file_tags = soup.find_all(node_type, {tag: re.compile(regex)})
    if not isinstance(filter, type(None)):
        file_tags = np.unique([filter(tag) for tag in file_tags], axis = 0).tolist()
    return file_tags

def crawl(urls, regex, filter_regex, session, label_filter = None, list_out = False):

    def _filter(tag):
        url = tag["href"]
        date_str = re.findall(filter_regex, url)[0].replace("/", "")
        return (date_str, url)

    if list_out:
        crawled_urls = list()
    else:
        crawled_urls = dict()

    for url in tqdm.tqdm(urls.values(), delay = 20):
        raw_paths = find_paths(url, regex, filter = _filter, session = session)
        if list_out:
            crawled_urls += raw_paths
        else:
            crawled_urls = {**crawled_urls, **{dt: new_url for dt, new_url in raw_paths if dt in label_filter}}

    return crawled_urls

def _download_urls(urls, folder, session, fps = None, parallel = 0, headers = None):

    if isinstance(parallel, bool):
        parallel = {True: -1, False: 0}[parallel]

    if isinstance(fps, type(None)):
        fps = [os.path.join(folder, os.path.split(url)[-1]) for url in urls]

    if parallel:
        backend = "loky"
        dler = partial(download_url, session = session, headers = headers)
        files = Parallel(n_jobs=parallel, backend = backend)(delayed(dler)(*x) for x in zip(urls, fps))
    else:
        files = list()
        for url, fp in zip(urls, fps):
            fp_out = download_url(url, fp, session = session, headers = headers)
            files.append(fp_out)

    return files

def update_urls(urls, dled_files, relevant_fps, checker_function = None):
    for fp in dled_files:
        if fp in relevant_fps:
            continue
        if isinstance(checker_function, type(None)):
            is_relevant = True
        else:
            is_relevant = checker_function(fp)
        l = [os.path.split(x)[-1] for x in urls]
        x = os.path.split(fp)[-1]
        url_idx = l.index(x) if x in l else None
        if os.path.isfile(fp):
            if not isinstance(url_idx, type(None)):
                _ = urls.pop(url_idx)
            if is_relevant:
                relevant_fps.append(fp)
    return urls, relevant_fps

def download_urls(urls, folder, session = None, fps = None, parallel = 0, headers = None, checker_function = None):
    
    relevant_fps = list()
    try_n = 0
    max_tries = 5
    url_fns = [os.path.split(x)[-1] for x in urls]
    try_again = True

    while try_again:

        try:
            _ = _download_urls(urls, folder, session, fps = fps, parallel = parallel, headers = headers)
        except NameError as e:
            log.info(f"--> {e}")
        finally:
            dled_files = glob.glob(os.path.join(folder, "*.nc"))
            dled_files = [x for x in dled_files if os.path.split(x)[-1] in url_fns]

        if not isinstance(checker_function, type(None)) and isinstance(fps, type(None)):
            try_n += 1
            urls, relevant_fps = update_urls(urls, dled_files, relevant_fps, checker_function)
            try_again = len(urls) > 0 and try_n < max_tries
        else:
            try_again = False

    if not isinstance(checker_function, type(None)):
        if len(urls) > 0:
            log.warning(f"--> Didn't succeed to download {len(urls)} files.")
        return relevant_fps
    else:
        return dled_files

def _download_url(url, fp, session = None, waitbar = True, headers = None):

    if os.path.isfile(fp):
        return fp

    folder, fn = os.path.split(fp)
    if not os.path.exists(folder):
        os.makedirs(folder)

    ext = os.path.splitext(fp)[-1]
    temp_fp = fp.replace(ext, "_temp")

    if isinstance(session, type(None)):
        session = requests.Session()

    file_object = session.get(url, stream = True, headers = headers)
    file_object.raise_for_status()

    if "Content-Length" in file_object.headers.keys():
        tot_size = int(file_object.headers["Content-Length"])
    else:
        tot_size = None

    if waitbar:
        wb = tqdm.tqdm(unit='Bytes', unit_scale=True, leave = False, 
                        total = tot_size, desc = fn)

    with open(temp_fp, 'wb') as z:
        for data in file_object.iter_content(chunk_size=1024):
            size = z.write(data)
            if waitbar:
                wb.update(size)

    os.rename(temp_fp, fp)

    return fp

def download_url(url, fp, session = None, waitbar = True, headers = None):
    max_tries = 10
    try_n = 0
    wait_sec = 15
    succes = False
    while try_n <= max_tries and not succes:
        if try_n > 0:
            waiter = int((try_n**1.2) * wait_sec)
            log.info(f"--> Trying to download {fp}, attempt {try_n+1} of {max_tries} in {waiter} seconds.")
            time.sleep(waiter)
        try:
            fp = _download_url(url, fp, session = session, waitbar = waitbar, headers = headers)
            succes = True
        except socket.timeout as e:
            log.info(f"--> Server connection timed out.")
        except HTTPError as e:
            log.info(f"--> Server error {e}.")
        else:
            ...
        finally:
            try_n += 1
    
    if not succes:
        raise NameError("Could not download {url} after {max_tries} attempts.")

    return fp

if __name__ == "__main__":

    import numpy as np

    folder = r"/Users/hmcoerver/On My Mac/crawl_test"
    _ = adjust_logger(True, folder, "INFO")

    hundreds = np.arange(101, 104)
    twohundreds = np.arange(200, 209)
    threehundreds = np.arange(300, 309)
    fourhundreds = np.arange(400, 430)
    fivehundreds1 = np.arange(500, 512)
    fivehundreds2 = np.arange(520, 528)

    all = np.concatenate([
                    # hundreds, 
                    twohundreds,
                    threehundreds,
                    fourhundreds,
                    fivehundreds1,
                    fivehundreds2])

    urls = [f"https://httpstat.us/{code}?test_{code}.txt" for code in all]

    folder1 = os.path.join(folder, "test1")
    if not os.path.isdir(folder1):
        os.makedirs(folder1)
    for fn in glob.glob(os.path.join(folder1, "*")):
        os.remove(fn)
    out1 = _download_urls(urls, folder1, requests.Session(), fps = None, parallel = 0, headers = None)