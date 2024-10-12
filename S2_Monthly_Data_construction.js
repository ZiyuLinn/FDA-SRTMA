//***********************************************************************
//    Function: Create monthly gap-filled features, surface reflectence, VI (NDVI, EVI, and EVI2)  for Sentinel-2 data
//       Input: Product information over flat surface (S2L2A has done the topographic and atmospheric corrections)
//   Reference: https://github.com/MarcYin/SIAC; https://eartharxiv.org/ps957/
//   Cloud removal: https://medium.com/google-earth/more-accurate-and-flexible-cloud-masking-for-sentinel-2-images-766897a9ba5f
// Acknowledge: Fengn Yin, UCL; ucfafy@ucl.ac.uk
//   Copyright: Shengbiao Wu @ HKU, 2020; Ziyu LIN @ HKU 2021
//***********************************************************************


// ----------------------------Notes--------------------------------------
// (REHHM is too time consuming, computational problem)
// gap-filling based on monthly median and then user SG filter for all bands and indies (9 window, 3 order)
// -----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Basic and default parameters
var s2cloud = ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY"),
    ecoregion = ee.FeatureCollection("users/linziyu/WesternGreatLake/northern_highland_level4");
ecoregion = ecoregion.geometry().transform('EPSG:4326', 10)//.aside(Map.addLayer)
var ecoregion_bound  = ecoregion.buffer(7000).bounds().aside(print)

Map.addLayer(ee.Feature(boundtiles.toList(6).get(2)))



var band_RF=ee.List([ 'B2', 'B3', 'B4','B5','B6','B7', 'B8', 'B8A','B11','B12']);
var band_VIs = ee.List(['NDVI','EVI','SAVI']);//,'DVI','NDRE','NDPI'
var all_bands = band_RF.cat(band_VIs)
var vizParams={bands:['B8','B4','B3'], min:0, max: 5000, gamma: [0.95, 1.1, 1]};
var vizParams_NDVI={bands:'NDVI', min:0, max: 1, gamma: [0.95]};
// var siac = require('users/marcyinfeng/Toa2Lai:S2_Toa2Lai');
// var siac = require('users/marcyinfeng/utils:SIAC');
var dem = ee.Image("USGS/3DEP/10m");
var palettes = require('users/gena/packages:palettes');
// var ltgee = require('users/emaprlab/public:Modules/LandTrendr.js'); // load the LandTrendr.js module

var lag = 20 // temporal resolution of the composites
var th = 0.5 // Define the percentage of amplitude for the estimation of the threshold
var YEAR = 2020 // year of processing
var start_date = YEAR+'-03-15'; // 5-10 is snow-free month , 3&11provide data for smoothing
var end_date   = YEAR+'-12-15';
print(start_date,end_date)

var PI = ee.Number(3.14159265359);
var MAX_SATELLITE_ZENITH = 6.0;
var MAX_DISTANCE = 1000000;
var UPPER_LEFT = 0;
var LOWER_LEFT = 1;
var LOWER_RIGHT = 2;
var UPPER_RIGHT = 3;
var MAX_CLOUD_PROBABILITY =50;  // threshold for cloud removal, set 50% here.
var Min_Snow_Threshold = 0.2 // minimum NDVI value for the reclassification of snow values


var SR_BAND_SCALE = 1e4
var CLOUD_FILTER = 50 //threshold for controlling cloud cover,The smaller the threshold, the more images will be excluded
var SNOW_FILTER = 50 //threshold for controlling snow cover
var CLD_PRB_THRESH = 10
var NIR_DRK_THRESH = 0.15
var B2_SNOW_THRESH = 0.1
var B11_SNOW_THRESH = 0.3
var CLD_PRJ_DIST = 2
var BUFFER = 20 // unit is meter
// Map.setZoom(14) // zoom level 14 = 10 meter resolution
// Map.setLocked(true, 14, 14)

//-----------------------------------Main function------------------------------------


//************************************************************************************************
// Step1: All available images preprocessing
// ************************************************************************************************
var S2_Level2A = ee.ImageCollection('COPERNICUS/S2_SR')


var n_tiles = boundtiles.size().getInfo()

var bound_list = boundtiles.toList(n_tiles+1).aside(print)
for(var i=2;i<3;i=i+1){

  if(ee.ImageCollection('projects/ee-project1-unmixing/assets/S2_ecoregion').filterMetadata('system:index','equals','tile_'+i).size().getInfo()){
    print('pass')
  continue
  }
  else{
  var site_bound = ee.Feature(bound_list.get(i)).geometry().buffer(20).bounds()
  var criteria4gapfill = ee.Filter.and(
    ee.Filter.bounds(site_bound), ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER),
      ee.Filter.lt('SNOW_ICE_PERCENTAGE', SNOW_FILTER),
    ee.Filter.date('2019-01-01', '2022-01-01'),
    ee.Filter.dayOfYear(75, 345)
    );
  var S2_collection_SR_4Gapfill = S2_Level2A
          .filter(criteria4gapfill)
          .map(addCloudMask)
          .map(function(img){ return img.clip(site_bound).reproject('EPSG:4326',null,10)})//.aside(print)
          .map(createDOYfield)
          .map(maskBadPixels)
          .map(add_cloud_bands).map(add_shadow_bands).map(add_cld_shdw_mask).map(apply_cld_shdw_mask)
          // this filter will only be conducted in winter (doy<150 or >250)
          .map(maskSnow)
          // .map(correct_SIAC)
          .map(BRDF_correction).map(PLC)
          .map(reflectence2float)
          // .map(addVIs)
          // .map(float2int)
          // .select(all_bands)
          .select(band_RF)
          .aside(print)

  var S2_collection_SR_4Gapfill_subtile = average_same_date(S2_collection_SR_4Gapfill).aside(print)



   var subtiles = site_bound.coveringGrid('EPSG:4326', 1000)//.aside(Map.centerObject,15)
  // exclude extreme value for each pixel
  // cloud removal according to time-series images
  var low = 5
  var up = 100-low
  var per_lst = subtiles.map(function(tile){
    var per_img = S2_collection_SR_4Gapfill_subtile.filterBounds(tile.geometry())
    .map(function(img){ return img.select('B2').clip(tile.geometry())}).reduce(ee.Reducer.percentile([low,up],['Plow','Pup']))
    return per_img
  })
  var per = ee.ImageCollection(per_lst).mosaic().clip(site_bound)
  // Map.addLayer(per)

  S2_collection_SR_4Gapfill_subtile = S2_collection_SR_4Gapfill_subtile
  .map(function(img){

  var cloud_mask = img.select('B2').gt(ee.Image(per).select(['.*_Plow'])).and(img.select('B2').lt(ee.Image(per).select(['.*Pup'])))
  .unmask(img.select('B2').mask())

  var cloud_mask_closing =  cloud_mask.eq(0)
  .focal_min(2).max(2).toInt()//.focalMin(1,'square',"pixels").focalMax(2,'square',"pixels")
  // .reduceResolution('mean') .reproject({'crs': img.select([0]).projection(), 'scale': 10}) // time-consuming

  // remove small cloud clusters inside the clear regions (usually are deciduous or bareland)
  // Remove isolated cloud clusters
  // connectedPixelCount is Zoom dependent, so visual result will vary
  var small_cloud = cloud_mask_closing.connectedPixelCount(16, true).lt(16)
  var cloud_mask2 = cloud_mask_closing.subtract(small_cloud).neq(1)


  // remove isolated small clear surface clusters outside the big clear regions
  var small_clear = cloud_mask2.connectedPixelCount(64, true).lt(64)
  var cloud_mask3 = cloud_mask2.subtract(small_clear).eq(1)

  var img_corr = img.updateMask(cloud_mask2.and(cloud_mask3))


  return img_corr//.addBands([cloud_mask.rename('m1'),cloud_mask_closing.rename('m2'),cloud_mask2.rename('m3'),cloud_mask3.rename('m4'),cloud_mask_closing.connectedPixelCount(64, true).rename('s1'),cloud_mask2.connectedPixelCount(64, true).rename('s2')])

  })



  //************************************************************************************************
  // Step2: filter out image for specific year and gap-filling by month
  // ************************************************************************************************
  // var S2_collection_SR_4Gapfill_smoothed= RMMEH(S2_collection_SR_4Gapfill, 7,all_bands)
  // Map.addLayer(S2_collection_SR_4Gapfill_smoothed)

  var S2_collection_SR = S2_collection_SR_4Gapfill_subtile//_smoothed
            .filterMetadata('Year','equals',YEAR)
            // .aside(print,'S2_collection_SR')


  // var S2_collection_SR_4Gapfill_ws= whittakerSmoothing(S2_collection_SR_4Gapfill)
  // Map.addLayer(S2_collection_SR_4Gapfill_ws)
  var month_lst = ee.List.sequence(3,12)//S2_collection_SR_4Gapfill_subtile.aggregate_array('Month').distinct()//.aside(print)
  var monthly_median_4gapfill =  ee.ImageCollection(month_lst.map(function(m){
    m = ee.Number(m).toInt()
    var median = S2_collection_SR_4Gapfill_subtile.filterMetadata('Month','equals',m)//.select(['.*RMMEH'])
    .median()
    // var bandnames = median.bandNames()
    return ee.Image(median).set('system:index',ee.String('month').cat(ee.String(m))).set('Month',m)//.reproject('EPSG:4326',null,10)
  }).flatten())
  // Bi-Monthly median
  var season_lst = month_lst.map(function(n){return ee.Number(n).divide(3).toInt()}).distinct()//.aside(print)
  var seasonal_median_4gapfill =  ee.ImageCollection(season_lst.map(function(s){
    s = ee.Number(s).toInt()
    var median = S2_collection_SR_4Gapfill_subtile.filter(ee.Filter.rangeContains('Month',s.multiply(3),s.add(1).multiply(3)))//.select(['.*RMMEH'])
    .median()
    // var bandnames = median.bandNames()
    return ee.Image(median).set('Season',s)//.reproject('EPSG:4326',null,10)
  }).flatten())
  var year_median = S2_collection_SR_4Gapfill.median()
  // Map.addLayer(monthly_median_4gapfill)


  //************************************************************************************************
  // Step3: gap fill nodata pixels using multi-years monthly avarage
  // ************************************************************************************************


  var S2_collection_SR_month_col = ee.ImageCollection(ee.List.sequence(3,12).map(function(m){
    m = ee.Number(m).toInt()
    var month_list = S2_collection_SR.filterMetadata('Month','equals',m).toList(30)
    var nan_img = year_median.mask(0)
    var month_list_nan = ee.List([nan_img.set('date',ee.Date.fromYMD(YEAR,m,1).millis()).set('system:time_start',ee.Date.fromYMD(YEAR,m,1).millis()).set('Month',m).set('Year',YEAR)
    ,nan_img.set('date',ee.Date.fromYMD(YEAR,m,15).millis()).set('system:time_start',ee.Date.fromYMD(YEAR,m,15).millis()).set('Month',m).set('Year',YEAR)
     ,nan_img.set('date',ee.Date.fromYMD(YEAR,m,28).millis()).set('system:time_start',ee.Date.fromYMD(YEAR,m,28).millis()).set('Month',m).set('Year',YEAR)])
    // make sure there are at least 3 image in one month
    var list = ee.Algorithms.If(month_list.size().eq(0),month_list_nan,
                ee.Algorithms.If(month_list.size().eq(1), month_list.cat(month_list_nan.slice(1,3)),
                 ee.Algorithms.If(month_list.size().eq(2), month_list.cat(month_list_nan.slice(2,3)),month_list)))
    return list
    }).flatten()).aside(print)

  var S2_col_gapfilled = S2_collection_SR_month_col.map(function(img){
    var month_median = monthly_median_4gapfill.filterMetadata('Month','equals',img.get('Month')).first()
    var season_median = seasonal_median_4gapfill.filterMetadata('Season','equals',ee.Number(img.get('Month')).divide(3).toInt()).first()
    // var new_bandnames = img.select(['.*_RMMEH']).bandNames().map(function(n){return ee.String(n).replace('_RMMEH','')})
    return img//.select(['.*_RMMEH'])
    .unmask(month_median,false)
    .unmask(season_median,false)
    .unmask(year_median,false)
  })

  //************************************************************************************************
  // Step3:  SG fitting by subtiles and output data
  // ************************************************************************************************
  // ============================= SG fitting  ====================================

  // Step 2: Set up Savitzky-Golay smoothing
  // Parameters: col, window, linear regression order
  // ** make sure have enough data within the window
  // (to solve the corr, valid observation should greater thant the linear regression order)
  // window = 7 is better, otherwise the spring will be strongly affected by summer
  // order 2,3 are similar
  var S2_col_sged =  SG_filter(S2_col_gapfilled, 7, 2)
  // Map.addLayer(S2_col_sged)

  var col = S2_col_sged.map(function(img){
    var int_img = img.multiply(10000).toUint16();
    return int_img.copyProperties(img,img.propertyNames().cat(['system:time_start']))})
    .map(createDOYfield)
    // .aside(print,'smoothed imagery')
  var month_lst = ee.List.sequence(4,11)//S2_col_sged.aggregate_array('Month').distinct()
  var monthly_median = month_lst.map(function(m){
    m = ee.Number(m).toInt()
    var month_col = col.filterMetadata('Month','equals',m).select('.*sg')
    // var median = ee.Algorithms.If(month_col.size().gt(0),month_col.median(),monthly_median_4gapfill.filterMetadata('Month','equals',m).first())
    var median = month_col.median()
    return ee.Image(median).set('system:index',ee.String('month').cat(ee.String(m)))
  }).flatten()

  var sg_QA = S2_collection_SR_4Gapfill.select('B2').count().rename('QA').toInt()//.eq(ee.List(month_lst).size())
  var img_monthly = ee.ImageCollection(monthly_median).toBands().addBands(sg_QA).toUint16();
  // Map.addLayer(img_monthly)
  site_bound.aside(Map.centerObject,14)
  // img_monthly = img_monthly.updateMask(img_monthly.mask().eq(1))
  Export.image.toAsset(img_monthly, 's2_monthly_tile'+i, 'projects/ee-project1-unmixing/assets/S2_ecoregion/tile_'+i,
  null, null, site_bound, 10, 'EPSG:4326',null,1e13)
  }

}

//********************************************************************************
// Function list
// Func 1: SIAC atmospheric correction
// function correct_SIAC(image)
// {
//   var cloud = s2cloud.filterMetadata('system:index', 'equals', image.get('system:index')).first()
//   // var atmcor = siac.inv_prosail(image);
//   var atmcor = siac.get_sur(image)
//   var boa=atmcor.select(band);
//   return boa;
// }
function average_same_date(col){
  var year_lst = col.aggregate_array('Year').distinct().map(function(y){
    var year_col = col.filterMetadata('Year','equals',y)
    var doy_lst = year_col.aggregate_array('DOY').distinct().map(function(d){
      var day_col = year_col.filterMetadata('DOY','equals',d)
      return ee.Image(day_col.first()).select([]).addBands(day_col.mean().rename(band_RF))
    })
    return doy_lst
  }).flatten()
  return ee.ImageCollection(year_lst)
}


// ===============Func 2: PLC topographic correction===============
function PLC(img){

// Extract image metadata about solar  and sensor postions
var SZ_rad = ee.Image.constant(ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE'))).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000));
var SA_rad = ee.Image.constant(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(3.14159265359).divide(180)).clip(img.geometry().buffer(10000));

var VZ_rad = ee.Image.constant(ee.Number(img.get('MEAN_INCIDENCE_ZENITH_ANGLE_B8'))).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000));
var VA_rad = ee.Image.constant(ee.Number(img.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8')).multiply(3.14159265359).divide(180)).clip(img.geometry().buffer(10000));

// Creat terrain layers
var slp = ee.Terrain.slope(dem).clip(img.geometry().buffer(10000));
var slp_rad = ee.Terrain.slope(dem).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000));
var asp_rad = ee.Terrain.aspect(dem).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000));


//var cosSolar=SZ_rad.cos().multiply(slp_rad.cos()).add(SZ_rad.sin().multiply(slp_rad.sin()).multiply(SA_rad.subtract(asp_rad)));
//var cosSensor=VZ_rad.cos().multiply(slp_rad.cos()).add(VZ_rad.sin().multiply(slp_rad.sin()).multiply(VA_rad.subtract(asp_rad)));
var cosSolar =SZ_rad.cos();
var cosSensor=VZ_rad.cos();
var cosSolart=SZ_rad.cos().multiply(ee.Image(1).subtract(slp_rad.tan().multiply(SZ_rad.tan()).multiply((SA_rad.subtract(asp_rad)).cos())));
var cosSensort=VZ_rad.cos().multiply(ee.Image(1).subtract(slp_rad.tan().multiply(VZ_rad.tan()).multiply((VA_rad.subtract(asp_rad)).cos())));

var up=(cosSolar.add(cosSensor)).multiply(cosSolart).multiply(cosSensort);
var down=(cosSolart.add(cosSensort)).multiply(cosSolar).multiply(cosSensor);
var factor_P=up.divide(down);

// var blue=img.select('B2').multiply(factor_P);
// var green=img.select('B3').multiply(factor_P);
// var red=img.select('B4').multiply(factor_P);
// var nir=img.select('B8').multiply(factor_P);

return img.multiply(factor_P).copyProperties(img,img.propertyNames().cat(['system:index','system:time_start']))

}

// ===============Func 3: BRDF correction for image collection===============
function BRDF_correction(image)
{
  var S2_SR_BRDF = applyBRDFMean_S2(image)//.multiply(0.0001); //BRDF L8
  return S2_SR_BRDF;
}

// Func 4: BRDF correction for single image, using the mean solar and sensor geometry
function applyBRDFMean_S2(image){
    var date = image.date();
    var sunAz = ee.Image(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(PI).divide(180));
    var sunZen = ee.Image(ee.Number(image.get('MEAN_SOLAR_ZENITH_ANGLE')).multiply(PI).divide(180));

    var viewAz = ee.Image(ee.Number(image.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8')).multiply(PI).divide(180));
    var viewZen = ee.Image(ee.Number(image.get('MEAN_INCIDENCE_ZENITH_ANGLE_B8')).multiply(PI).divide(180));

    var kval_Normal = _kvol(ee.Image(0).multiply(PI.divide(180)), ee.Image(45).multiply(PI.divide(180)), ee.Image(0).multiply(PI.divide(180)), ee.Image(0).multiply(PI.divide(180)));
    var kval_Sensor = _kvol(sunAz, sunZen, viewAz, viewZen);
    var result = _applyS2(image, kval_Normal,kval_Sensor);

    return result//.addBands(sunZen).addBands(sunAz).addBands(viewZen).addBands(viewAz);
}


// calculate the normal solar zenith usign the approach of Claverie et al. 2018
function _sunZen_Normal(footprint){
  var lat_centor = ee.Number(ee.List((footprint.get(0))).get(1)).add(ee.List(footprint.get(2)).get(1)).divide(2)
  var k = ee.List([6.15e-11,-1.95e-09,-9.48e-7,2.4e-5,0.0119,-0.127,31])
  var sunZen_n = lat_centor.pow(6).multiply(value(k,0))
  .add(lat_centor.pow(5).multiply(value(k,1)))
  .add(lat_centor.pow(4).multiply(value(k,2)))
  .add(lat_centor.pow(3).multiply(value(k,3)))
  .add(lat_centor.pow(2).multiply(value(k,4)))
  .add(lat_centor.pow(1).multiply(value(k,5))).add(value(k,6))
  return ee.Image(sunZen_n)
  //var sunZen_n = lat_centor.multiply()
}

// get solar zenith angle from satellite (L8&S2) images
// refs: https://www.sciencedirect.com/science/article/abs/pii/S0960148104003404
function getsunAngles(date, footprint){
  var jdp = date.getFraction('year');
  var seconds_in_hour = 3600;
  var  hourGMT = ee.Number(date.getRelative('second', 'day')).divide(seconds_in_hour);
  var latRad = ee.Image.pixelLonLat().select('latitude').multiply(PI.divide(180));
  var longDeg = ee.Image.pixelLonLat().select('longitude');

  // Julian day proportion in radians
  var jdpr = jdp.multiply(PI).multiply(2);

  var a = ee.List([0.000075, 0.001868, 0.032077, 0.014615, 0.040849]);
  var meanSolarTime = longDeg.divide(15.0).add(ee.Number(hourGMT)); // Local solar time, Lon/15+hourGMT
  var localSolarDiff1 = value(a, 0)
          .add(value(a, 1).multiply(jdpr.cos()))
          .subtract(value(a, 2).multiply(jdpr.sin()))
          .subtract(value(a, 3).multiply(jdpr.multiply(2).cos()))
          .subtract(value(a, 4).multiply(jdpr.multiply(2).sin()));
  //E
  var localSolarDiff2 = localSolarDiff1.multiply(12 * 60);
  //E*720
  var localSolarDiff = localSolarDiff2.divide(PI);  // E*720/PI
  var trueSolarTime = meanSolarTime
          .add(localSolarDiff.divide(60))
          .subtract(12.0);
  // Lon/15+hourGMT+E*720/PI/60-12
  // Hour as an angle;
  var ah = trueSolarTime.multiply(ee.Number(MAX_SATELLITE_ZENITH * 2).multiply(PI.divide(180))) ;
  var b = ee.List([0.006918, 0.399912, 0.070257, 0.006758, 0.000907, 0.002697, 0.001480]);
  //(Lon/15+hourGMT+E*720/PI/60-12)*15*PI/180
  var delta = value(b, 0)
        .subtract(value(b, 1).multiply(jdpr.cos()))
        .add(value(b, 2).multiply(jdpr.sin()))
        .subtract(value(b, 3).multiply(jdpr.multiply(2).cos()))
        .add(value(b, 4).multiply(jdpr.multiply(2).sin()))
        .subtract(value(b, 5).multiply(jdpr.multiply(3).cos()))
        .add(value(b, 6).multiply(jdpr.multiply(3).sin()));

  var cosSunZen = latRad.sin().multiply(delta.sin())
        .add(latRad.cos().multiply(ah.cos()).multiply(delta.cos()));
  var sunZen = cosSunZen.acos();

  // sun azimuth from south, turning west
  var sinSunAzSW = ah.sin().multiply(delta.cos()).divide(sunZen.sin());
  sinSunAzSW = sinSunAzSW.clamp(-1.0, 1.0);

  var cosSunAzSW = (latRad.cos().multiply(-1).multiply(delta.sin())
                    .add(latRad.sin().multiply(delta.cos()).multiply(ah.cos())))
                    .divide(sunZen.sin());
  var sunAzSW = sinSunAzSW.asin();

  sunAzSW = where(cosSunAzSW.lte(0), sunAzSW.multiply(-1).add(PI), sunAzSW);
  sunAzSW = where(cosSunAzSW.gt(0).and(sinSunAzSW.lte(0)), sunAzSW.add(PI.multiply(2)), sunAzSW);

  var sunAz = sunAzSW.add(PI);

   // # Keep within [0, 2pi] range
    sunAz = where(sunAz.gt(PI.multiply(2)), sunAz.subtract(PI.multiply(2)), sunAz);

  var footprint_polygon = ee.Geometry.Polygon(footprint);
  sunAz = sunAz.clip(footprint_polygon);
  sunAz = sunAz.rename(['sunAz']);
  sunZen = sunZen.clip(footprint_polygon).rename(['sunZen']);

  return [sunAz, sunZen];
}

// get azimuth angle from satelltie (L8&S2) image
function azimuth(footprint){
    function x(point){return ee.Number(ee.List(point).get(0))}
    function  y(point){return ee.Number(ee.List(point).get(1))}

    var upperCenter = line_from_coords(footprint, UPPER_LEFT, UPPER_RIGHT).centroid().coordinates();
    var lowerCenter = line_from_coords(footprint, LOWER_LEFT, LOWER_RIGHT).centroid().coordinates();
    var slope = ((y(lowerCenter)).subtract(y(upperCenter))).divide((x(lowerCenter)).subtract(x(upperCenter)));
    var slopePerp = ee.Number(-1).divide(slope);
    var azimuthLeft = ee.Image(PI.divide(2).subtract((slopePerp).atan()));
    return azimuthLeft.rename(['viewAz']);
  }

function zenith(footprint){
    var leftLine = line_from_coords(footprint, UPPER_LEFT, LOWER_LEFT);
    var rightLine = line_from_coords(footprint, UPPER_RIGHT, LOWER_RIGHT);
    var leftDistance = ee.FeatureCollection(leftLine).distance(MAX_DISTANCE);
    var rightDistance = ee.FeatureCollection(rightLine).distance(MAX_DISTANCE);
    var viewZenith = rightDistance.multiply(ee.Number(MAX_SATELLITE_ZENITH * 2))
          .divide(rightDistance.add(leftDistance))
          .subtract(ee.Number(MAX_SATELLITE_ZENITH))
          .clip(ee.Geometry.Polygon(footprint))
          .rename(['viewZen']);

    //print(footprint)
    //print(leftDistance)
    return viewZenith.multiply(PI.divide(180));
}

// BRDF correction for whole image
function _applyS2(image, kval_normal, kval_sensor){
      var f_iso = 0;
      var f_geo = 0;
      var f_vol = 0;

    //don't have parameters for red edge. Use the parameters from nearest red or nir band instead
			var blue = _correct_band(image, 'B2', kval_normal, kval_sensor, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372);
			var green = _correct_band(image, 'B3', kval_normal, kval_sensor, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580);
			var red = _correct_band(image, 'B4', kval_normal, kval_sensor, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574);
      var nir = _correct_band(image, 'B8', kval_normal, kval_sensor, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535);
      var rededge1 = _correct_band(image, 'B5', kval_normal, kval_sensor, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574);
      var rededge2 = _correct_band(image, 'B6', kval_normal, kval_sensor, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574);
      var rededge3 = _correct_band(image, 'B7', kval_normal, kval_sensor, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535);
      var rededge4 = _correct_band(image, 'B8A', kval_normal, kval_sensor, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535);
      var swir1 = _correct_band(image, 'B11', kval_normal, kval_sensor, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154);
      var swir2 = _correct_band(image, 'B12', kval_normal, kval_sensor, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639);
			return image.select([]).addBands([blue, green, red, rededge1,rededge2,rededge3,nir,rededge4,swir1,swir2])
}

// Band specific BRDF correction
function _correct_band(image, band_name, kval_normal, kval_sensor, f_iso, f_geo, f_vol){
	//"""fiso + fvol * kvol + fgeo * kgeo"""
	var iso = ee.Image(f_iso);
	var geo = ee.Image(f_geo);
	var vol = ee.Image(f_vol);
  var kvol_n = kval_normal[0]
  var kgeo_n = kval_normal[1]
  var kvol_s = kval_sensor[0]
  var kgeo_s = kval_sensor[1]
	var pred0 = vol.multiply(kvol_n).add(geo.multiply(kgeo_n)).add(iso).rename(['pred0']);
	var pred = vol.multiply(kvol_s).add(geo.multiply(kgeo_s)).add(iso).rename(['pred']);
	var cfac = pred0.divide(pred).rename(['cfac']);
	var corr = image.select(band_name).multiply(cfac).rename([band_name]);
	return corr;
}

// Kernel function
function _kvol(sunAz, sunZen, viewAz, viewZen){
	//"""Calculate kvol kernel.
	//From Lucht et al. 2000
	//Phase angle = cos(solar zenith) cos(view zenith) + sin(solar zenith) sin(view zenith) cos(relative azimuth)"""

  var relative_azimuth = sunAz.subtract(viewAz).rename(['relAz']);

  var cos_xi = viewZen.cos()
    .multiply(sunZen.cos())
    .add(viewZen.sin()
    .multiply(sunZen.sin())
    .multiply(relative_azimuth.cos()));

  var xi = cos_xi.acos();

  var kvol = (((xi.multiply(-1).add(Math.PI/2))
    .multiply(cos_xi)
    .add(xi.sin())).divide(sunZen.cos().add(viewZen.cos())))
    .subtract(Math.PI/4).rename(['kvol']);

  var h_b = 2;
  var b_r = 1;

  var theta_ip = (sunZen.tan().multiply(b_r)).atan();
  var theta_vp = (viewZen.tan().multiply(b_r)).atan();

  var D = (theta_ip.tan().pow(2)
    .add(theta_vp.tan().pow(2))
    .subtract(theta_ip.tan()
    .multiply(theta_vp.tan())
    .multiply(relative_azimuth.cos())
    .multiply(2))).sqrt();

  var cos_t1 = (D.pow(2)
    .add((theta_ip.tan()
    .multiply(theta_vp.tan())
    .multiply(relative_azimuth.sin())).pow(2)))
    .sqrt()
    .multiply(h_b)
    .divide(theta_ip.cos().pow(-1).add(theta_vp.cos().pow(-1)))

  var cos_t2 = cos_t1.where(cos_t1.lt(-1),-1)
  var cos_t = cos_t2.where(cos_t2.gt(1),1)

  var t = cos_t.acos()

  var O = (t.subtract(t.sin().multiply(cos_t)))
    .multiply(theta_ip.cos().pow(-1).add(theta_vp.cos().pow(-1)))
    .divide(Math.PI)

  var cos_xip = theta_ip.cos()
    .multiply(theta_vp.cos())
    .add(theta_ip.sin()
    .multiply(theta_vp.sin())
    .multiply(relative_azimuth.cos()));

  var kgeo = O.subtract(theta_ip.cos().pow(-1))
    .subtract(theta_vp.cos().pow(-1))
    .add((cos_xip.add(1)).multiply(viewZen.cos().pow(-1)).multiply(sunZen.cos().pow(-1)).multiply(0.5));

	return [kvol,kgeo]
}



// ===============Func 5: adds a doy propterty to the image.===============
function createDOYfield(image) {
  var date_ = ee.Date(image.get('system:time_start'))
  // var first_date = ee.Date.fromYMD(date_.get('year').subtract(1), 12, 31)
  var year = date_.get('year')
  var first_date = ee.Date.fromYMD(ee.Number(year).subtract(1), 12, 31)
  var doy = date_.difference(first_date,'day').int()
  var month = doy.divide(30).add(1).int()
  return image.set('DOY',doy.toInt()).set('Month',month.toInt()).set('Year',year)
.set('date', date_.millis()).set('oid',image.get('system:index'))
}


// ===============Func 6: adds a doy band to the image.===============
function addDataBands(image) {
  // return image.addBands(image.metadata('system:time_start').divide(1e18).rename('time'));
  var date_ = ee.Date(image.get('system:time_start'))
  var first_date = ee.Date.fromYMD(date_.get('year').subtract(1), 12, 31)
  var doy = date_.difference(first_date,'day').int()
  return image.addBands(ee.Image(doy).rename('DOY').cast({'DOY':'double'}))
}

// ===============Func 7: RMMEH smooth approach=============
// ref: Jin, Z., & Xu, B. (2013). A novel compound smoother—RMMEH to reconstruct MODIS NDVI time series. IEEE Geoscience and Remote Sensing Letters, 10(4), 942-946.
// max(avrage of two neigbor, median of the whole window, center)
// @ collection -- can contains several bands that need to be smoothed
// @ window -- must be singal int
// @ bandNames -- list of bands that need to be smoothed
function RMMEH(collection, window, bandNames){
  // collection: colection to smooth - must be single band
  // window: smooth window. Time that will be added forward and back
  // window = 5 will give a [t-2, t-1, t, t+1, t+1] window for the medians smoother
  // window = 3 will give a [t-1, t, t+1] window for the medians smoother
 var half_window = (window - 1)/2;
 bandNames = ee.List(bandNames).sort()
 var new_bandNames = bandNames.map(function(n){return ee.String(n).cat('_RMMEH')}).sort()
 collection=collection.sort('DOY').select(bandNames);
// var collection_list=collection.toList(collection.size());
 var doy_list = collection.aggregate_array('DOY').distinct().sort()
 // slice -- strat(inclusive), end(exclusive)
// var doy_first = doy_list.slice(half_window,half_window+1).get(0);
// var doy_last = doy_list.slice(-half_window,-half_window+1).get(0);

  // create a list collection based on DOYS
  // This approach aims to reduce the computation memory
  // var col_list = ee.List(doy_list.iterate(function(n,list){
  //   return ee.List(list).add(ee.ImageCollection(collection.filterMetadata('DOY','equals',ee.Number(n))))},
  //   ee.List([])))
  //   .aside(print,'col_list')
  var col_list = ee.ImageCollection(doy_list.map(function(n){
    var col_doy = collection.filterMetadata('DOY','equals',ee.Number(n))
    return col_doy.median().set('DOY',n).set('system:time_start',ee.Date('2020-01-01').advance(n,'day'))
     })).sort('DOY').toList(366)
    .aside(print,'col_list')

  var col_list_extend = col_list.slice(0,half_window+1).reverse().cat(col_list)
  .cat(col_list.slice(doy_list.size().subtract(half_window),doy_list.size().add(1)))
  print(col_list.get(0))


  // Run the medians over the window, and return the smoothed image version
  // var image_first=ee.ImageCollection([ee.Image(collection_list.get(ee.Number(0))),ee.Image(collection_list.get(ee.Number(1)))]);
  // var image_last=ee.ImageCollection([ee.Image(collection_list.get(ee.Number(-2))),ee.Image(collection_list.get(ee.Number(-1)))]);
  // var image_first=collection.filterMetadata('DOY','not_greater_than',doy_first).sort('DOY').toList(20);
  // var image_last=collection.filterMetadata('DOY','not_less_than',doy_last).sort('DOY').toList(20);
  // collection_list = image_last.cat(collection_list).cat(image_first.reverse())
  // strat(inclusive), end(inclusive)
  var runLength = ee.List.sequence(0, doy_list.size().subtract(1));
  var median_series = runLength.map(function(n) {
    var num_center = ee.Number(n).add(half_window)
    var image_center = ee.Image(col_list_extend.get(num_center));
    // var doy_center= ee.Number(doy_list.get(n))
    // var image_center = ee.ImageCollection(col_list_extend.get(num_center))
    // var image_center_median = image_center.median();

    var image_left = col_list_extend.slice(num_center.subtract(half_window),num_center)
    var image_right = col_list_extend.slice(num_center.add(1),num_center.add(half_window+1))
    var col_range = ee.ImageCollection(col_list_extend.slice(num_center.subtract(half_window),num_center.add(half_window+1)))
    var image_median=col_range.median();
    var image_average=ee.ImageCollection([ee.ImageCollection(image_left).mean(),ee.ImageCollection(image_right).mean()]).mean();

    // var doy_center = ee.Number(image_center.get('DOY'))
    //  (inclusive)
    // var col = collection.filter(ee.Filter.rangeContains('DOY',doy_center.subtract(half_window),doy_center.add(half_window)))
    // List.slice: Start value must be less than or equal to end value.
    // var col= ee.ImageCollection(col_list_extend.slice(num_center.subtract(half_window-1),num_center.add(half_window+1))
    // .iterate(function(l,first){ return ee.ImageCollection(first).merge(l)}, col_list_extend.get(num_center.subtract(half_window))))
    // // var image_left = collection.filter(ee.Filter.rangeContains('DOY',doy_center.subtract(half_window),doy_center.subtract(1)))
    // // var image_right =collection.filter(ee.Filter.rangeContains('DOY',doy_center.add(1),doy_center.add(half_window)))
    // var image_median=col.median();
    // var image_average=col.mean();


    // extend window to generate the median value for gap-filling
    // var image4gapfill  = ee.ImageCollection(collection_list.slice(ee.Algorithms.If(num_center.subtract(window).lt(0),0,num_center.subtract(window)),
    // ee.Algorithms.If(num_center.add(window).gt(collection.size()),collection.size(),num_center.add(window)))).median()

    // if image are missing during this period (for several years), them the value is probably contaiminated (snowy or rainy month)
    var valid_mask = col_range.select([bandNames.get(0)]).count()

    var image_all=ee.ImageCollection([image_center,image_median,image_average]).median()//.mask(valid_mask)//.unmask(image4gapfill)//.toFloat();


    return image_center.addBands([image_all.rename(new_bandNames),valid_mask.divide(10).rename('valid_counts')])
    // .copyProperties(image_center.sort('system:time_start').first(),image_center.sort('system:time_start').first().propertyNames().cat(['system:index','system:time_start']));

    // return image_center.select([]).addBands(image_all.rename(new_bandNames)).copyProperties(image_center,image_center.propertyNames().cat(['system:index','system:time_start']));

  });
  return  ee.ImageCollection(median_series); //image_first.merge(ee.ImageCollection(median_series)).merge(image_last).sort("DOY");
}

//===============Func 8: Sentinel S2 cloud mask===========================
// The masks for the 10m bands sometimes do not exclude bad data at
// scene edges, so we apply masks from the 20m and 60m bands as well.
function maskBadPixels(s2_img) {
  return s2_img.updateMask(
      s2_img.select('B8A').mask()
      // .updateMask(s2_img.select('B9').mask())
      );
}

function add_cloud_bands(img){
    // Get s2cloudless image, subset the probability band.
    var cld_prb = ee.Image(img).select('probability')

    // Condition s2cloudless by the probability threshold value.
    var is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds')

    // Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([is_cloud]))
}


function add_shadow_bands(img){
    // Identify water pixels from the SCL band.
    // var not_water = img.select('SCL').neq(6)

    // for level 1 data, do not have SCL layer, use NDWI for water dectection
    var nir = img.select("B8");
    var green = img.select("B3");
    var not_water = img.expression(
  "(B3 - B8)/(B3 + B8)",
  { "B8": nir, "B3": green}).lt(0.5);

    // Identify dark NIR pixels that are not water (potential cloud shadow pixels).

    var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')

    // Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    // Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    var cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject({'crs': img.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform'))

    // Identify the intersection of dark pixels with cloud shadow projection.
    var shadows = cld_proj.multiply(dark_pixels).rename('shadows')

    // Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))
}

function add_cld_shdw_mask(img){
    // Add cloud component bands.
    var img_cloud = add_cloud_bands(img)

    // Add cloud shadow component bands.
    var img_cloud_shadow = add_shadow_bands(img_cloud)

    // Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    var is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)

    // Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    // 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    var is_cld_shdw = (is_cld_shdw.focal_min(2).focal_max(BUFFER*2/20)
        .reproject({'crs': img.select([0]).projection(), 'scale': 20})
        .rename('cloudmask'))

    // Add the final cloud-shadow mask to the image.
    return img_cloud_shadow.addBands(is_cld_shdw)
}

function apply_cld_shdw_mask(img){
    // Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
    var not_cld_shdw = ee.Image(img).select('cloudmask').not()

    // Subset reflectance bands and update their masks, return the result.
    return img.select('B.*').updateMask(not_cld_shdw)
}
function addCloudMask(image){
  var col = s2cloud.filterMetadata('system:index', 'equals', image.get('system:index'))
  // if no cloud mask, cloud probability ==100, label as useless image
  var cloudmask = ee.Algorithms.If(col.size().eq(0),ee.Image(100).rename('probability'),col.first())

  return image.addBands(cloudmask)
          }

// Function to remove cloud and snow pixels
// only apply on the image in winter
function maskSnow(image){
  // Normalised Difference Snow Index
  // Source: https://earth.esa.int/web/sentinel/technical-guides/sentinel-2-msi/level-2a/algorithm
  // NDSI above 0.2 are usually snow,0.2~0.42 is partly snowy, <0.2 is non-snow, water sometimes has high NDSI (~0)
  var snow = image.normalizedDifference(['B3','B11']).lt(0.2)
  .and(image.select(['B2']).lt(B2_SNOW_THRESH*SR_BAND_SCALE))
  .and(image.select(['B11']).lt(B11_SNOW_THRESH*SR_BAND_SCALE)).focal_min(3).focal_max(3)
  var doy = ee.Number(image.get('DOY'))

  return ee.Algorithms.If(doy.gt(270).or(doy.lt(150)),image.updateMask(snow),image);
}


// ============Func 9: add Vegetation indies=========================
function reflectence2float(img){
  return ee.Image(img).divide(10000).toFloat().copyProperties(img,img.propertyNames().cat(['system:index','system:time_start']))
}
function float2int(img){
  return ee.Image(img).multiply(10000).toInt16().copyProperties(img,img.propertyNames().cat(['system:index','system:time_start']))
}
function addVIs(img){
  var image = ee.Image(img)
  var NDVI = image.normalizedDifference(['B5', 'B4']).rename('NDVI').float()
  // var NDWI = image.normalizedDifference(['B3', 'B8']).rename('NDWI').float()
  var EVI = image.expression(
    '2.5 * ((nir - red) / (nir + 6 * red - 7.5 * blue + 1))', {
      'nir' : image.select('B8'),
      'red' : image.select('B4'),
      'blue' : image.select('B2')
    }).rename('EVI').float();
  var SAVI = image.expression(
    '((nir-red) / (nir+red+0.428)) * (1.428)', {
      'nir' : image.select('B8'),
      'red' : image.select('B4')
    }).rename('SAVI').float();
  // var DVI = image.expression(
  //   '(edgeTwo-red)', {
  //     'edgeTwo' : image.select('B6'),
  //     'red' : image.select('B4')
  //   }).rename('DVI').float();
  // // normalized difference red-edge
  // var NDRE =  image.expression(
  //   '((nir-edgeTwo) / (nir+edgeTwo))', {
  //     'nir' : image.select('B8'),
  //     'edgeTwo' : image.select('B6')
  //   }).rename('NDRE').float();
    // define a new red-SWIR reflectance as the weighted sum of the red and SWIR reflectance
    // in order to minimize the difference be- tween soil and snow
  // var NDPI = image.expression(
  //   '(nir-(0.74*red+0.26*swir))/(nir+(0.74*red+0.26*swir))',{
  //     'nir' : image.select('B8'),
  //     'red' : image.select('B4'),
  //     'swir' : image.select('B11')
  //   }).rename('NDPI').float();
  return image.addBands([NDVI,EVI,SAVI])
}

function VIs_mask(img){
var mask = ee.Image(img).select(band_VIs).reduce('min').gt(0);
var img_ = img.updateMask(mask)
var isWinter = ee.Number(img.get('DOY')).mod(365)

var img_ = ee.Algorithms.If(isWinter.gt(300).or(isWinter.lt(100)),img_.updateMask(ee.Image(img).select('NDPI').lt(0.8)),img_)
return ee.Image(img_).copyProperties(img)
}


//====================Func 10:  Whittaker smoother:====================
// Kong, D., Zhang, Y., Gu, X., & Wang, D. (2019). A robust method for reconstructing global MODIS EVI time series on the Google Earth Engine. *ISPRS Journal of Photogrammetry and Remote Sensing*, *155*(May), 13–24.
// https://doi.org/10.1016/j.isprsjprs.2019.06.014
// Note: gee_whittaker is designed for MODIS VIs. Sentinel-2 might be less satisfactory.

// helper function to convert qa bit image to flag
function extractBits(image, start, end, newName) {
    // Compute the bits we need to extract.
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += Math.pow(2, i);
    }
    // Return a single band image of the extracted QA bits, giving the band
    // a new name.
    return image.select([0], [newName])
                  .bitwiseAnd(pattern)
                  .rightShift(start);
}

// function to get a Difference mattrix of specified order
// on the input matrix. takes matrix and order as parameters
function getDifferenceMatrix(inputMatrix, order){
    var rowCount = ee.Number(inputMatrix.length().get([0]));
    var left = inputMatrix.slice(0,0,rowCount.subtract(1));
    var right = inputMatrix.slice(0,1,rowCount);
    if (order > 1 ){
        return getDifferenceMatrix(left.subtract(right), order-1)}
    return left.subtract(right);
};

// unpacks an array image into images and bands
// takes an array image, list of image IDs and list of band names as arguments
function unpack(arrayImage, imageIds, bands){

    function iter(item, icoll){

        function innerIter(innerItem, innerList){
            return ee.List(innerList).add(ee.String(item).cat("_").cat(ee.String(innerItem)))}

        var temp = bands.iterate(innerIter, ee.List([]));
        return ee.ImageCollection(icoll)
            .merge(ee.ImageCollection(ee.Image(arrayImage).select(temp,bands).set("id",item)))}

    var imgcoll  = ee.ImageCollection(imageIds.iterate(iter, ee.ImageCollection([])));
    return imgcoll
}

//Function to compute the inverse log ratio of a regression results to
// transform back to percent units
function inverseLogRatio(image) {
  var bands = image.bandNames();
  var t = image.get("system:time_start");
  var ilrImage = ee.Image(100).divide(ee.Image(1).add(image.exp())).rename(bands);
  return ilrImage.copyProperties(image,image.propertyNames().cat(['system:index','system:time_start']))
}

function whittakerSmoothing(imageCollection,isCompositional, lambda){
  // quick configs to set defaults
  if (isCompositional === undefined || isCompositional !==true) isCompositional = false;
  if (lambda === undefined ) lambda = 10;

  // procedure start
  var ic = imageCollection.map(function(image){
     var props = image.toDictionary();
    return image.toFloat().set(props);
  }).aside(print);

  var dimension = ic.size();
  var identity_mat = ee.Array.identity(dimension);
  var difference_mat = getDifferenceMatrix(identity_mat,3);
  var difference_mat_transpose = difference_mat.transpose();
  var lamda_difference_mat = difference_mat_transpose.multiply(lambda);
  var res_mat = lamda_difference_mat.matrixMultiply(difference_mat);
  var hat_matrix = res_mat.add(identity_mat);


  // backing up original data
  var original = ic;

  // get original image properties
  var properties = ee.List(ic.iterate(function(image, list){
    return ee.List(list).add(image.toDictionary());
  },[]));

  var time = ee.List(ic.iterate(function(image, list){
    return ee.List(list).add(image.get("system:time_start"));
  },[]));

  // if data is compositional
  // calculate the logratio of an image between 0 and 100. First
  // clamps between delta and 100-delta, where delta is a small positive value.
  if (isCompositional){
    ic = ic.map(function(image){
      var props = image.toDictionary();
      var delta = 0.001;
      var bands = image.bandNames();
      image = image.clamp(delta,100-delta);
      image = (ee.Image.constant(100).subtract(image)).divide(image).log().rename(bands);
      return image.set(props);
    });
  }

  var arrayImage = original.toArray();
  var coeffimage = ee.Image(hat_matrix);
  var smoothImage = coeffimage.matrixSolve(arrayImage);

  var idlist = ee.List(ic.iterate(function(image, list){
    return ee.List(list).add(image.get("system:index"));
  },[])).aside(print,'idlist');
  var bandlist = ee.Image(ic.first()).bandNames();

  var flatImage = smoothImage.arrayFlatten([idlist,bandlist]);
  var smoothCollection = ee.ImageCollection(unpack(flatImage, idlist, bandlist)).aside(print,'wsmoothed');

  if (isCompositional){
    smoothCollection = smoothCollection.map(inverseLogRatio);
  }
  // get new band names by adding suffix fitted
  var newBandNames = bandlist.map(function(band){return ee.String(band).cat("_fitted")});
  // rename the bands in smoothened images
  smoothCollection = smoothCollection.map(function(image){return ee.Image(image).rename(newBandNames)});

  // a really dumb way to loose the google earth engine generated ID so that the two
  // images can be combined for the chart
  var dumbimg = arrayImage.arrayFlatten([idlist,bandlist]);
  var dumbcoll = ee.ImageCollection(unpack(dumbimg,idlist, bandlist));
  var outCollection = dumbcoll.combine(smoothCollection);

  var outCollectionProp = outCollection.iterate(function(image,list){
    var id = idlist.get(ee.List(list).size());
    var t = time.get(ee.List(list).size());
    return ee.List(list).add(image.set("system:index",id).set('system:time_start',t));
  },[]);

  // var outCollectionProp = outCollection.iterate(function(image,list){
  //   image = image.set("system:time_start",time.get(ee.List(list).size())).set()
  //   return ee.List(list).add(image);
  // },[]);


  var residue_sq = smoothImage.subtract(arrayImage).pow(ee.Image(2)).divide(dimension);
  // var rmse_array = residue_sq.arrayReduce(ee.Reducer.sum(),[0]).pow(ee.Image(1/2));

  // var rmseImage = rmse_array.arrayFlatten([["rmse"],bandlist]);

  // return [ee.ImageCollection.fromImages(outCollectionProp), rmseImage];
   return ee.ImageCollection.fromImages(outCollectionProp);
}

//============= Function: logistic function fitting for SG filter ================

// sg filter transcribed from:
// http://www2.geog.ucl.ac.uk/~plewis/eoldas/plot_directive/savitzky_golay.py
// Author: Guido Lemoine - EC JRC, 2017-02-23  transcribed from python library
//                                 2018-06-17  modified to compare to matrixSolve solution
//                                 2018-06-19  adapted to run on image collections
// Part 1: Set up a S1 collection that covers the aoi consistently.
// Best to select a limited area that is covered by a single orbit.
// If the selection has holes, the smoothing will throw errors (for now)
function SG_filter(col, window,order){
  var bands = ee.Image(col.first()).bandNames()
  var half_window = (window - 1)/2
  // Define the axes of variation in the collection array.
  var imageAxis = 0;
  var bandAxis = 1;
  var nbands = 1;


  // Set polynomial order
  if(order == 3){
  var coeffFlattener = [['constant', 'x', 'x2', 'x3']]
  var indepSelectors = ['constant', 't', 't2', 't3']
  var col_res = col.map(function(img) {
  var ddiff = ee.Date(img.get('date')).difference(ee.Date(start_date), 'hour')
    return img.addBands(ee.Image(1).toFloat().rename('constant'))
      .addBands(ee.Image(ddiff).toFloat().rename('t'))
      .addBands(ee.Image(ddiff).pow(ee.Image(2)).toFloat().rename('t2'))
      .addBands(ee.Image(ddiff).pow(ee.Image(3)).toFloat().rename('t3'))
      .updateMask(img.mask().reduce('min'))
  }).sort('date')
  }
  else if(order = 2){
  var coeffFlattener = [['constant', 'x', 'x2']]
  var indepSelectors = ['constant', 't', 't2']
   var col_res = col.map(function(img) {
  var ddiff = ee.Date(img.get('date')).difference(ee.Date(start_date), 'hour')
  return img.addBands(ee.Image(1).toFloat().rename('constant'))
      .addBands(ee.Image(ddiff).toFloat().rename('t'))
      .addBands(ee.Image(ddiff).pow(ee.Image(2)).toFloat().rename('t2'))
      .updateMask(img.mask().reduce('min'))
  }).sort('date')
  }


  // For the remainder, use s1res as a list of images
  // set time stamp for regression and sorting

  var col_list = col_res.toList(col_res.size())
  var runLength = ee.List.sequence(0, col_res.size().subtract(half_window+1))//.aside(print,'runLength')
  // Run the SG solver over the series, and return the smoothed image version
  var sg_series = runLength.map(function(i) {
    var ref = ee.Image(col_list.get(ee.Number(i).add(half_window)))
    var h = ref.get('date')
    // Convert the collection to an array.
    var img_list = bands.map(function(b){
      var array_ = col_res.select(ee.List([b]).cat(indepSelectors)).toArray();

     // Get a slice corresponding to the window_size of the SG smoother
      var subarray = array_.arraySlice(imageAxis, ee.Number(i).int(), ee.Number(i).add(window).int())
      var predictors = subarray.arraySlice(bandAxis, nbands, nbands + order + 1)
      var response = subarray.arraySlice(bandAxis, 0, 1); // VV
      var coeff = predictors.matrixSolve(response)

      coeff = coeff.arrayProject([0]).arrayFlatten(coeffFlattener)

      var sged = coeff.multiply(ref.select(indepSelectors)).reduce(ee.Reducer.sum()).rename([ee.String(b).cat('_sg')])
      return ref.select([b]).addBands(sged)
    })
    if(bands.size().gt(1)){
    var img_sged = ee.Image(img_list.slice(1,bands.size()).iterate(function(i,f){ return ee.Image(f).addBands(i)},img_list.get(0)))
    }else{ var img_sged = ee.Image(img_list.get(0))}
    return img_sged.copyProperties(ref).set('system:index',ee.String(h))
  })

  return ee.ImageCollection(sg_series)
}
