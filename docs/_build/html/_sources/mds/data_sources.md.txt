To run the ETLook model, two types of spatial variables are required, temporal and static data. **Each of these variables can be collected from whichever source you wish to use**, as long as you make sure the units are correct, the data is stored as a GeoTIFF (1 band per file, 1 file for each variable and date), the files all have the same no-data-value and they all have the same projection and resolution.

**For your convenience, the pyWAPOR package has a function that can collect all this data from selected sources** and make sure the data is stored in the correct format and folder structure.

#### Temporal ET_Look Data (composites)
| Variable | Unit | Selected Sources |
| ------ | ------ | ------ |
| Normalized Difference Vegetation Index (NDVI) | — | MOD13, MYD13, PROBA-V|
| Albedo | — | MCD43, PROBA-V|
| Precipitation | mm/day | CHIRPS |
| Air Pressure at sea level | kPa | MERRA-2, GEOS-5 |
| Specific Humidity | kg/kg | MERRA-2, GEOS-5 |
| Air Temperature | °C | MERRA-2, GEOS-5 |
| Windspeed | m/s | MERRA-2, GEOS-5 |
| Solar Radiation | W/m2  | MERRA-2 |
| Soil Saturation | — | from pywapor.se_root() |

#### Temporal SE_Root Data (instantaneous)
| Variable | Unit | Selected Sources |
| ------ | ------ | ------ |
| Normalized Difference Vegetation Index (NDVI) | — | MOD13, MYD13, PROBA-V |
| Air Pressure at sea level | kPa | MERRA-2, GEOS-5 |
| Air Pressure at surface level  | kPa | MERRA-2, GEOS-5 |
| Specific Humidity | kg/kg | MERRA-2, GEOS-5 |
| Air Temperature  | °C | MERRA-2, GEOS-5 |
| Windspeed | m/s | MERRA-2, GEOS-5 |
| Total Precipitable Water Vapour  | mm | MERRA-2, GEOS-5 |
| Land Surface Temperature (LST) | K | MOD11, MYD11 |

#### Static Data
| Variable | Unit | Selected Sources |
| ------ | ------ | ------ |
| Landcover | — | WaPOR, GlobCover |
| Digital Elevation Model | m.a.s.l | SRTM |
| Air Temperature (yearly amplitude) | K | GLDAS |
| Latitude | DD | from NDVI |
| Longitude | DD | from NDVI |
| Slope | ° | from Elevation |
| Slope Aspect | ° | from Elevation |
| Bulk Stomatal Resistance | s/m | from Landcover |
| Landmask | — | from Landcover |
| Maximum Light Use Efficiency | gr/MJ | from Landcover |
| Maximum Obstacle Height | m | from Landcover |

#### Sources
| Source | Temporal Availability | Temporal Resolution |Spatial Resolution | Used For |
| ------ | ------ | ------ | ------ | ------ |
|[MOD13](https://lpdaac.usgs.gov/products/mod13q1v006/) | 2000-02-18 - ongoing | 16-Daily |250m|NDVI|
|[MYD13](https://lpdaac.usgs.gov/products/myd13q1v006/) | 2002-07-04 - ongoing | 16-Daily |250m|NDVI|
|[MCD43](https://lpdaac.usgs.gov/products/mcd43a3v006/)|2000-02-16 - ongoing|Daily|500m|Albedo|
|[MOD11](https://lpdaac.usgs.gov/products/mod11a1v006/) | 2000-02-24 - ongoing | Daily | 1000m | LST |
|[MYD11](https://lpdaac.usgs.gov/products/myd11a1v006/)| 2002-07-04 - ongoing | Daily | 1000m | LST |
|[PROBAV](https://www.vito-eodata.be/collectioncatalogue/srv/eng/catalog.search#/metadata/urn:ogc:def:EOP:VITO:PROBAV_S5-TOC_100M_V001)|2014-03-11 - ongoing|5-Daily|100m|NDVI, Albedo|
| [GEOS5](https://geos5.org) | 2017-12-01 - ongoing | 3-Hourly |0.3125°×0.25° | Meteo |
| [MERRA2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/) | 1980-01-01 - ongoing | Hourly | 0.625°×0.5° | Meteo | 
| [CHIRPS](https://www.chc.ucsb.edu/data/chirps) |  1981-01-01 - ongoing | Daily | 0.05° | Precipitation |
| [WAPOR](https://wapor.apps.fao.org/catalog/WAPOR_2/1/L1_LCC_A) | 2009 - 2020 | Yearly |250m | Landcover |
| [GLOBCOVER](http://due.esrin.esa.int/page_globcover.php) | 2009 | Single| 250m | Landcover |
| [SRTM](https://srtm.csi.cgiar.org) | 2009 | Single | 90m | DEM |