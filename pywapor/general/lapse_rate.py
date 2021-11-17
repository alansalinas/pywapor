import pywapor.general.processing_functions as PF
import numpy as np
import pywapor

def lapse_rate_temperature(tair_file, dem_file):

    # import matplotlib.pyplot as plt

    # tair_file = r"/Volumes/Data/pre_et_look_NEW/RAW/MERRA/Air_Temperature/daily_MERRA2/t2m_MERRA_K_daily_2019.07.06.tif"
    # tair_file = r"/Volumes/Data/pre_et_look_NEW/RAW/GEOS/Air_Temperature/daily/t2m_GEOS_K_daily_2019.07.06.tif"

    # def _plot_array(array):
    #     plt.clf()
    #     plt.imshow(array)
    #     plt.title(array.shape)
    #     plt.colorbar()
    #     plt.show()

    ## 1
    ds_t_down = PF.reproject_dataset_example(tair_file, dem_file, 2)
    tempe = PF.open_as_array(ds_t_down)
    # _plot_array(tempe - 273.15)

    ## 2
    dem_down = PF.open_as_array(dem_file)
    # _plot_array(dem_down)

    ## 3
    ds_dem_up = PF.reproject_dataset_example(dem_file, tair_file, 4)
    
    dem_up = PF.open_as_array(ds_dem_up)
    dem_up[np.isnan(dem_up)] = 0.
    ds_dem_up = PF.Save_as_MEM(dem_up, ds_dem_up.GetGeoTransform(), PF.Get_epsg(ds_dem_up))
    
    ds_dem_up_down = PF.reproject_dataset_example(ds_dem_up, dem_file, 2)
    dem_up_ave = PF.open_as_array(ds_dem_up_down)
    # _plot_array(dem_up_ave)

    ## Correct wrong values
    dem_down[dem_down <= 0] = 0
    dem_up_ave[dem_up_ave <= 0] = 0

    tdown = pywapor.et_look_v2.meteo.disaggregate_air_temperature(tempe, dem_down, dem_up_ave, lapse = pywapor.et_look_v2.constants.lapse)
    # _plot_array(tdown)

    fh = tair_file.replace("_K_", "_C_")
    geo, projection = PF.get_geoinfo(dem_file)[:2]

    PF.Save_as_tiff(fh, tdown, geo, projection)

    # test_tair = r"/Volumes/Data/pre_et_look_ORIGINAL/ETLook_input_MODIS/20190706/tair_24_20190706.tif"
    # tempe_test = PF.open_as_array(test_tair)
    # _plot_array(tempe_test - tdown)

    return fh