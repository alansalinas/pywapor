from sentinelsat import SentinelAPI, read_geojson, geojson_to_wkt
from datetime import date

api = SentinelAPI('broodj3ham', 'nuqso2-wicxIv-wyxcix', 'https://apihub.copernicus.eu/apihub')

footprint = geojson_to_wkt(read_geojson(r"/Users/hmcoerver/Desktop/inst3_2d_asm_Nx.geojson"))
products = api.query(footprint,
                     date=(date(2021, 8, 29), date(2021, 12, 29)),
                     platformname='Sentinel-2',
                     cloudcoverpercentage=(0, 40))