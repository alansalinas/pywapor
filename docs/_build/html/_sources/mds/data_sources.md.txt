To run the ETLook model, two types of spatial variables are required, temporal and static data. Each of these variables can be collected from whichever source you wish to use, as long as you make sure the units are correct and the projection consistent. .

For your convenience, the pyWaPOR package has a function that can collect the required data from selected sources and make sure the data is stored in the correct format and folder structure.

If you are interested in using data from sources not listed below, please have a look at the "sideloading" notebook in the User Guide.

<p><details open>
<summary><b>Data sources</b></summary>
<br>
<table class = "docutils align-default">
    <thead>
        <tr class="row-odd" style="text-align:center">
            <th class="head">Source</td>
            <th class="head">Temporal Availability</td>
            <th class="head">Temporal Resolution</td>
            <th class="head">Spatial Resolution</td>
            <th class="head">Used For</td>
        </tr>
    <thead>
    <tr class="row-odd">
        <td><a href="https://lpdaac.usgs.gov/products/mod13q1v006/">MOD13</a></td>
        <td>2000-02-18 - ongoing</td>
        <td>16-Daily</td>
        <td>250m</td>
        <td>NDVI</td>
    </tr>
    <tr class="row-odd">
        <td><a href="https://lpdaac.usgs.gov/products/myd13q1v006/">MYD13</a></td>
        <td>2002-07-04 - ongoing</td>
        <td>16-Daily</td>
        <td>250m</td>
        <td>NDVI</td>
    </tr>
    <tr class="row-odd">
        <td><a href="https://lpdaac.usgs.gov/products/mcd43a3v006/">MCD43</a></td>
        <td>2000-02-16 - ongoing</td>
        <td>Daily</td>
        <td>500m</td>
        <td>Albedo</td>
    </tr>
    <tr class="row-odd">
        <td><a href="https://lpdaac.usgs.gov/products/mod11a1v006/">MOD11</a></td>
        <td>2000-02-24 - ongoing</td>
        <td>Daily</td>
        <td>1000m</td>
        <td>LST</td>
    </tr>
    <tr class="row-odd">
        <td><a href="https://lpdaac.usgs.gov/products/myd11a1v006/">MYD11</a></td>
        <td>2002-07-04 - ongoing</td>
        <td>Daily</td>
        <td>1000m</td>
        <td>LST</td>
    </tr>
    <tr class="row-odd">
        <td><a href="https://www.vito-eodata.be/collectioncatalogue/srv/eng/catalog.search#/metadata/urn:ogc:def:EOP:VITO:PROBAV_S5-TOC_100M_V001">PROBAV</a></td>
        <td>2014-03-11 - ongoing</td>
        <td>5-Daily</td>
        <td>100m</td>
        <td>NDVI, Albedo</td>
    </tr>
    <tr class="row-odd">
        <td><a href="https://geos5.org">GEOS5</td>
        <td>2017-12-01 - ongoing</a></td>
        <td>3-Hourly</td>
        <td>0.3125°×0.25°</td>
        <td>Meteo</td>
    </tr>
    <tr class="row-odd">
        <td><a href="https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/">MERRA2</a></td>
        <td>1980-01-01 - ongoing</td>
        <td>Hourly</td>
        <td>0.625°×0.5°</td>
        <td>Meteo</td>
    </tr>
    <tr class="row-odd">
        <td><a href="https://www.chc.ucsb.edu/data/chirps">CHIRPS</a></td>
        <td>1981-01-01 - ongoing</td>
        <td>Daily</td>
        <td>0.05°</td>
        <td>Precipitation</td>
    </tr>
    <tr class="row-odd">
        <td><a href="https://wapor.apps.fao.org/catalog/WAPOR_2/1/L1_LCC_A">WAPOR</a></td>
        <td>2009 - 2020</td>
        <td>Yearly</td>
        <td>250m</td>
        <td>Landcover</td>
    </tr>
    <tr class="row-odd">
        <td><a href="http://due.esrin.esa.int/page_globcover.php">GLOBCOVER</a></td>
        <td>2009</td>
        <td>Single</td>
        <td>250m</td>
        <td>Landcover</td>
    </tr>
    <tr class="row-odd">
        <td><a href="https://srtm.csi.cgiar.org">SRTM</a></td>
        <td>2009</td>
        <td>Single</td>
        <td>90m</td>
        <td>DEM</td>
    </tr>
</table>

</details></p>

<p><details>
<summary><b>Temporal data (composites)</b></summary>
<br>
<table class = "docutils align-default">
   <thead>
      <tr class="row-odd" style="text-align:center">
         <th class="head">Variable</th>
         <th class="head">Unit</th>
         <th class="head">Selected Sources</th>
      </tr>
   </thead>
   <tbody>
      <tr class="row-odd">
         <td>Normalized difference vegetation index</td>
         <td>-</td>
         <td>MOD13, MYD13, PROBAV</td>
      </tr>
      <tr class="row-even">
         <td>Albedo</td>
         <td>-</td>
         <td>MCD43, PROBAV</td>
      </tr>
      <tr class="row-odd">
         <td>Precipitation</td>
         <td>mm/day</td>
         <td>CHIRPS</td>
      </tr>
      <tr class="row-even">
         <td>Air pressure at sea level</td>
         <td>mbar</td>
         <td>MERRA2, GEOS5</td>
      </tr>
      <tr class="row-odd">
         <td>Specific humidity</td>
         <td>kg/kg</td>
         <td>MERRA2, GEOS5</td>
      </tr>
      <tr class="row-even">
         <td>Air temperature</td>
         <td>°C</td>
         <td>MERRA2, GEOS5</td>
      </tr>
      <tr class="row-odd">
         <td>Windspeed</td>
         <td>m/s</td>
         <td>MERRA2, GEOS5</td>
      </tr>
      <tr class="row-even">
         <td>Solar radiation</td>
         <td>W/m<sup>2</sup></td>
         <td>MERRA2</td>
      </tr>
      <tr class="row-odd">
         <td>Soil saturation</td>
         <td>-</td>
         <td>pywapor.se_root()</td>
      </tr>
   </tbody>
</table>

</details></p>

<p><details>
<summary><b>Temporal data (instantaneous)</b></summary>
<br>
<table class = "docutils align-default">
    <thead>    
        <tr class="row-odd" style="text-align:center">
            <th class="head">Variable</td>
            <th class="head">Unit</td>
            <th class="head">Selected Sources</td>
        </tr>
    </thead>
    <tr class="row-odd">
        <td>Normalized difference vegetation index</td>
        <td>—</td>
        <td>MOD13, MYD13, PROBAV</td>
    </tr>
    <tr class="row-even">
        <td>Air pressure at sea level</td>
        <td>kPa</td>
        <td>MERRA2, GEOS5</td>
    </tr>
    <tr class="row-odd">
        <td>Air pressure at surface level</td>
        <td>kPa</td>
        <td>MERRA2, GEOS5</td>
    </tr>
    <tr class="row-even">
        <td>Specific humidity</td>
        <td>kg/kg</td>
        <td>MERRA2, GEOS5</td>
    </tr>
    <tr class="row-odd">
        <td>Air temperature</td>
        <td>°C</td>
        <td>MERRA2, GEOS5</td>
    </tr>
    <tr class="row-even">
        <td>Windspeed</td>
        <td>m/s</td>
        <td>MERRA2, GEOS5</td>
    </tr>
    <tr class="row-odd">
        <td>Total precipitable water vapour</td>
        <td>mm</td>
        <td>MERRA2, GEOS5</td>
    </tr>
    <tr class="row-even">
        <td>Land surface temperature</td>
        <td>K</td>
        <td>MOD11, MYD11</td>
    </tr>
</table>

</details></p>

<p><details>
<summary><b>Static data</b></summary>
<br>
<table class = "docutils align-default">
    <thead>
        <tr class="row-odd" style="text-align:center">
            <th class="head">Variable</td>
            <th class="head">Unit</td>
            <th class="head">Selected Sources</td>
        </tr>
    </thead>
    <tr class="row-odd">
        <td>Landcover</td>
        <td>—</td>
        <td>WaPOR, GlobCover</td>
    </tr>
    <tr class="row-even">
        <td>Elevation</td>
        <td>m.a.s.l</td>
        <td>SRTM</td>
    </tr>
    <tr class="row-odd">
        <td>Air temperature (yearly amplitude)</td>
        <td>K</td>
        <td>GLDAS</td>
    </tr>
    <tr class="row-even">
        <td>Latitude</td>
        <td>DD</td>
        <td>from NDVI</td>
    </tr>
    <tr class="row-odd">
        <td>Longitude</td>
        <td>DD</td>
        <td>from NDVI</td>
    </tr>
    <tr class="row-even">
        <td>Slope</td>
        <td>°</td>
        <td>from Elevation</td>
    </tr>
    <tr class="row-odd">
        <td>Slope aspect</td>
        <td>°</td>
        <td>from Elevation</td>
    </tr>
    <tr class="row-even">
        <td>Bulk Stomatal Resistance</td>
        <td>s/m</td>
        <td>from Landcover</td>
    </tr>
    <tr class="row-odd">
        <td>Landmask</td>
        <td>—</td>
        <td>from Landcover</td>
    </tr>
    <tr class="row-even">
        <td>Maximum Light Use Efficiency</td>
        <td>gr/MJ</td>
        <td>from Landcover</td>
    </tr>
    <tr class="row-odd">
        <td>Maximum Obstacle Height</td>
        <td>m</td>
        <td>from Landcover</td>
    </tr>
</table>

</details></p>