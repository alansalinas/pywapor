import os 
import urllib
from bs4 import BeautifulSoup
import urllib
from joblib import Parallel, delayed
from functools import partial
import re
import numpy as np
import tqdm
import requests
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

def download_urls(urls, folder, session, fps = None, parallel = 0, headers = None):

    if isinstance(parallel, bool):
        parallel = {True: -1, False: 0}[parallel]

    if isinstance(fps, type(None)):
        fps = [os.path.join(folder, os.path.split(url)[-1]) for url in urls]

    dler = partial(download_url, session = session, headers = headers)

    if parallel:
        backend = "loky"
        files = Parallel(n_jobs=parallel, backend = backend)(delayed(dler)(*x) for x in zip(urls, fps))
    else:
        files = [dler(*x) for x in zip(urls, fps)]

    return files

def download_url(url, fp, session = None, waitbar = True, headers = None):

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