import json
from landsatxplore.api import API

# Initialize a new API instance and get an access key
api = API("Bert1989", "macmus-0cyqKe-qyjbyx")

# Search for Landsat TM scenes
scenes = api.search(
    dataset='landsat_ot_c2_l2',
    latitude=50.85,
    longitude=-4.35,
    start_date='2020-01-01',
    end_date='2022-03-01',
    max_cloud_cover=10
)

print(f"{len(scenes)} scenes found.")

# Process the result
for scene in scenes:
    print(scene['acquisition_date'])
    # Write scene footprints to disk
    # fname = f"{scene['landsat_product_id']}.geojson"
    # with open(fname, "w") as f:
    #     json.dump(scene['spatialCoverage'], f)

api.logout()

# from landsatxplore.earthexplorer import EarthExplorer

# ee = EarthExplorer("Bert1989", "macmus-0cyqKe-qyjbyx")

# ee.download(scene["entity_id"], output_dir='/Users/hmcoerver/Downloads/merra2')
