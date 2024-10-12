//***********************************************************************
// Function: monthly gap-filled features for Sentinel-1 data
//  Copyright: Ziyu LIN @ HKU 2022
//***********************************************************************



// ========================My updates ===============================
// add terrainCorrection, change dem to ee.Image("USGS/3DEP/10m") 2022.01.03
// Smooth results 2022.01.03


var ecoregion = ee.FeatureCollection("users/linziyu/WesternGreatLake/northern_highland_level4");
var aoi =  ee.FeatureCollection("users/linziyu/WesternGreatLake/northern_highland_level4")
ecoregion = aoi.geometry().transform('EPSG:4326', 10)//.aside(Map.addLayer)
var ecoregion_bound  = ecoregion.buffer(7000).bounds().aside(print).aside(Map.addLayer)
var boundtiles = ecoregion_bound.coveringGrid('EPSG:4326',10000)
var start_date = '2018-01-01'
var end_date = '2023-01-01'

// get the data from Sentinel-1 collection, for our area and period of interest
// Imagery in the Earth Engine 'COPERNICUS/S1_GRD' Sentinel-1 ImageCollection is consists of
// Level-1 Ground Range Detected (GRD) scenes processed to backscatter coefficient (σ°) in decibels (dB)
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterMetadata('instrumentMode', 'equals', 'IW')
   .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
  .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV', 'VH']))
  .filter(ee.Filter.dayOfYear(75,345))
  .filterMetadata('resolution_meters','equals',10)
  .filterBounds(aoi)
  .filterDate(start_date, end_date)
  .sort('system:time_start')




// Make sure to convert from dB(current unit) to "natural" values (power scale) before doing any calculations
// https://hyp3-docs.asf.alaska.edu/guides/introduction_to_sar/#power-scale
var n_tiles = boundtiles.size().getInfo()
var bound_list = boundtiles.toList(n_tiles+1).aside(print)
for(var i=0;i<142;i=i+1){
  if(ee.ImageCollection('projects/ee-project1-unmixing/assets/S1_ecoregion_v2').filterMetadata('system:index','equals','s1_tile'+i).size().getInfo()){
    print('pass')
  continue
  }
  else{
  var site_bound = ee.Feature(bound_list.get(i)).geometry().buffer(20).bounds().aside(print).aside(Map.centerObject)
  // Map.addLayer(site_bound)
  // Map.centerObject(site_bound,12)

  // Step 1: filter images and clip; exclue invalid data and extreme outliners
  var s1_col = s1.filterBounds(site_bound)
  .map(function(img){ return img.reproject('EPSG:4326',null,10).clip(site_bound)})
  // .aside(print)//.aside(Map.addLayer,{min:-15,max:0,bands:['VV']})
  .map(setDstamp)//.aside(print)
  .map(terrainCorrection)//.aside(Map.addLayer,{min:-15,max:0,bands:['VV']})//.aside(print)
  .map(angleNormalization)//.aside(Map.addLayer,{min:-15,max:0,bands:['VV']})//.aside(print)
  .map(toNatural)
  .map(RefinedLee)
  .map(RVI)
  .map(function(img){return img.updateMask(img.select('mRFDI').gt(0.1).and(img.gt(0)).and(img.lt(1)))})
  // .aside(print);

  // comparing filtered vs non-filtered
  // s1_col.aside(Map.addLayer,{min:-15,max:0})
  // s1_col.map(GammaMap).aside(Map.addLayer,{min:-15,max:0})
  // s1_col.map(RefinedLee).aside(Map.addLayer,{min:-15,max:0})


  // Step 2: merge image on the same time interval and gap-fill by monthly and yearly data
  var month_lst = s1_col.aggregate_array('month').distinct()
  var s1_month_median = ee.ImageCollection(month_lst
  .map(function(m){ return s1_col.filterMetadata('month','equals',m).median().set('month',m)}))
  var year_median = s1_col2.median()
  // Step 3: mosaic images in the same date and gap-fill using the monthly median
  var time_interval = 'week'
  var time_lst = s1_col.aggregate_array(time_interval).distinct()
  var s1_mosaic = ee.ImageCollection(time_lst.map(function (n) {
  var same_col = s1_col2.filterDate('2020-01-01', '2020-12-31').filterMetadata(time_interval, 'equals', n);
  var m = s1_col2.filterMetadata(time_interval, 'equals', n).first().get('month');
  var month_median = s1_month_median.filterMetadata('month', 'equals', m).first();
  // make sure there is image for every time stamp
  //better to use monthly median for gap-filling for weekly images for the lack of weekly valid data
  var week_median = ee.Algorithms.If(same_col.size().gt(0), same_col.median(), s1_col2.filterMetadata('month', 'equals', m).median());
  var dstamp = ee.Date('2019-12-28').advance(n, time_interval);
  // when doing unmask to mosaic, make sure to generate a new image instead of using the old one (image.select)
  // because the system:footprint properties will be used to constraint the image boundary
  // and make sure to set the unmask(sameFootprint=false)
  return ee.Image(week_median).unmask(month_median, false).unmask(year_median, false)  // .updateMask(year_mean.mask().reduce('min'))
.set(time_interval, n).set('month', m).set('system:time_start', dstamp.millis()).set('date', dstamp.millis()).set('oid', ee.String(time_interval).cat(n));
}))
// .aside(print,'mosaic')
// Map.addLayer(s1_mosaic)

  // Run the SG solver over the series, and return the smoothed image version
  var order =2
  var window_size = 7 //Big window will flattern the result if there is only a short time period peak
  // make sure the data are adequate for the 2-order regression
  // Solve the tile error: matrixSolve cannot solve underdetermined system.
  // If there are fewer equations than variables, then the system is called underdetermined
  var s1_QA = s1_col.select('VV').count().rename('QA')//.aside(Map.addLayer,{max:6,min:8})

  // Run the SG solver over the series, and return the smoothed image version
  var sg_series = SG_filter(s1_mosaic,time_interval, window_size, order)
  // Map.addLayer(sg_series)

  var s1_monthly_lst = ee.ImageCollection(sg_series).aggregate_array('month').distinct()
  .map(function(m){
  m = ee.Number(m).toInt()
  var mean = ee.ImageCollection(sg_series).filterMetadata('month','equals',m).median()

  return ee.Image(mean).set('system:index',ee.String('month').cat(ee.String(m)))
}).flatten()//.aside(print)
var s1_monthly_img = ee.ImageCollection(s1_monthly_lst).toBands().multiply(10000).toInt16().select('.*_sg').addBands(s1_QA).aside(print)


Export.image.toAsset(s1_monthly_img, 'img_s1sged_monthly_tile'+i, 'projects/ee-project1-unmixing/assets/S1_ecoregion_v2/s1_tile'+i,
null, null, site_bound, 10, 'EPSG:4326',null,1e13)
}
}



// ============================= Functions ===========================
function angleNormalization(img){
    //  the normalised sentinel-1 global backscatter model, mapping earth's land surface with c-band microwaves
    //  https://doi.org/10.3390/rs13152856,  backscatter intensity normalization based on Cosine model (norm to 38 degree)
    var normalized_ind = ee.Image(40).multiply(Math.PI / 180).cos().divide(img.select('angle').multiply(Math.PI / 180).cos()).pow(2)
    //exclude angle band
    // should be apply in the natural/linear format
    var img_norm = ee.Image(toDB(ee.Image(toNatural(ee.Image(img).select(['..']))).multiply(normalized_ind)));
    // var vmask = img_norm.select('VV').subtract(img_norm.select('VH')).expression('b(0) >0? 1:0');
    //  masking low dB returns (water/and maybe little bareland),
    // masking high dB returns (high buildings/snow)
    // var veg_mask = img_norm.expression('((VV<-15)|(VH<-22.9))?0:((VV>-2)?0:1)',{'VV':img_norm.select('VV'),'VH':img_norm.select('VH') });
    var snow_mask = img_norm.expression('((VV>-2)?0:1)',{'VV':img_norm.select('VV'),'VH':img_norm.select('VH') });
    var output = ee.Algorithms.If(ee.Number(img.get('month')).lt(5).or(ee.Number(img.get('month')).gt(9)),img.select([]).addBands(img_norm),
    img.select([]).addBands(img_norm).updateMask(snow_mask))
    return img.select([]).addBands(output)
  }





// Radar Vegetation Index Code for Dual Polarimetric Script
// https://custom-scripts.sentinel-hub.com/sentinel-1/radar_vegetation_index_code_dual_polarimetric/#
function RVI(img){
  /*
Dual-Pol Radar Vegetation index for Sentinel-1

This code is based on:
[1]RVI: Nasirzadehdizaji, Rouhollah, et al. "Sensitivity Analysis of Multi-Temporal Sentinel-1 SAR Parameters to Crop Height and Canopy Coverage." Applied Sciences 9.4 (2019): 655.
https://doi.org/10.3390/app9040655
[2]RFVI:Nicolau, A. P., Flores-Anderson, A., Griffin, R., Herndon, K., & Meyer, F. J. (2021). Assessing SAR C-band data to effectively distinguish modified land uses in a heavily disturbed Amazon forest. International Journal of Applied Earth Observation and Geoinformation, 94, 102214.
https://doi.org/10. 1016/j.jag.2020.102214
[3]DpRVI:Dual polarimetric radar vegetation index for crop growth monitoring using sentinel-1 SAR data. https://doi.org/10.1016/j.rse.2020.111954

*/

// formula:  (4xVH)/(VV+VH)
// must be in linear units ??? -> conversion formula: 10^((valueINdb)/10) DONE


  img = ee.Image(img)
  // var img_linear = ee.Image(toNatural(img.select(['..'])))
  //depolarization within the vegetation
  // var value = img.expression('((4*VH)/(VV+VH))',{'VV':img.select('VV'),'VH': img.select('VH')})
  // .rename('RVI');
  // var value2 = img.expression('((VV/(VV+VH))**0.5*(4*VH)/(VV+VH))',{'VV':img.select('VV'),'VH': img.select('VH')})
  // .rename('mRVI');
  var value3 = img.expression('((VV-VH)/(VV+VH))',{'VV':img.select('VV'),'VH': img.select('VH')})
  .rename('mRFDI');
  // var value4 = img_linear.expression('((VV/(VV+VH))**0.5*(4*VH)/(VV+VH))',{'VV':img_linear.select('VV'),'VH': img_linear.select('VH')})
  // .rename('mRVI_natural');
  // var value5 = img_linear.expression('((VV-VH)/(VV+VH))',{'VV':img_linear.select('VV'),'VH': img_linear.select('VH')})
  // .rename('mRFDI_natural');
  var dpRVI = DpRVI(img)

  // var entropy = img.select(['..']).toInt32().entropy(ee.Kernel.circle({radius:1.5,units:'pixels'}))
  // return img.addBands([value,value2,value3,img_linear,value4,value5,dpRVI])
  return img.addBands([value3])
}

/*----------------------------------------------------------------------------------------------
                      3) Generating Dual-pol descriptors and the clusters
----------------------------------------------------------------------------------------------*/


function DpRVI(image) {
    var window_size = 1 // 3*3
    var C11_mean = image.select('VV')//.expression( '10 ** (VV / 10)', {'VV': image.select('VV')})
                  .reduceNeighborhood({
                    reducer: ee.Reducer.mean(),
                    kernel: ee.Kernel.square(window_size)
                    });
    var C22_mean = image.select('VH')//.expression( '10 ** (VH / 10)', {'VH': image.select('VH')})
                  .reduceNeighborhood({
                    reducer: ee.Reducer.mean(),
                    kernel: ee.Kernel.square(window_size)
                    });
    // var C11_mean = image.select('VV')
    // var C22_mean = image.select('VH')

    var span = C11_mean.add(C22_mean);
    var ratio = C22_mean.divide(C11_mean);
    // var vmask = C11_mean.subtract(C22_mean);
    // vmask = vmask.expression('b(0) >0? 1:0');

    var m = (C11_mean.subtract(C22_mean).abs()).divide(span);
    var d_dpol = m.multiply(m).subtract(1).multiply(-1);
    var theta_c = ((C11_mean.subtract(C22_mean).abs()).multiply(span).multiply(m))
                        .divide((C11_mean.multiply(C22_mean)).add(span.pow(2).multiply(m.pow(2))))
                        .atan();
    // theta_c = theta_c.multiply(180).divide(Math.PI);

    var p1 = C11_mean.divide(span);
    var p2 = C22_mean.divide(span);
    var cnst = ee.Number(2);
    var Hp1 = p1.multiply(p1.log10()).divide(cnst.log10()).multiply(-1);
    var Hp2 = p2.multiply(p2.log10()).divide(cnst.log10()).multiply(-1);
    var H = Hp1.add(Hp2);
    var q = ratio;
    var DpRVIc_n = q.multiply(q.add(ee.Number(3)));
    var DpRVIc_d = (q.add(ee.Number(1))).multiply(q.add(ee.Number(1)));
    var DpRVIc = DpRVIc_n.divide(DpRVIc_d);
    // var C11_mean_db = C11_mean.log10().multiply(10);//Linear to dB conversion
    // var C11_rc = C11_mean_db.expression('b(0)<-17?0:1'); // masking low dB returns (water)



    //Masked values

    var out_raster = m.addBands([//theta_c, m,H
                                DpRVIc]);//.updateMask(vmask)//.updateMask(C11_rc);

    out_raster = out_raster.select(
        out_raster.bandNames(), // old names
        ['mc','DpRVIc']
        );
    return out_raster.set('system:time_start', image.get('system:time_start'));
};



// set time stamp for regression and sorting
function setDstamp(img){
var dstamp = ee.Date(img.get('system:time_start'))
var year = dstamp.get('year')
var month = dstamp.get('month')
var week = dstamp.get('week')
 var first_date = ee.Date.fromYMD(year.subtract(1), 12, 31)
  var doy = dstamp.difference(first_date,'day').int()
return img.set('date', dstamp).set('DOY',doy.toInt()).set('month',month.toInt()).set('week',week.toInt()).set('oid',ee.String(img.get('system:index')))
}

// Solve
function getLocalFit(i,array) {
  // Get a slice corresponding to the window_size of the SG smoother
  var subarray = array.arraySlice(imageAxis, ee.Number(i).int(), ee.Number(i).add(window_size).int())
  var predictors = subarray.arraySlice(bandAxis, nbands, nbands + order + 1)
  var response = subarray.arraySlice(bandAxis, 0, 1); // VV
  var coeff = predictors.matrixSolve(response)

  coeff = coeff.arrayProject([0]).arrayFlatten(coeffFlattener)
  return coeff
}

// Functions to convert from/to dB
function toNatural(img) {
  return img.select([]).addBands(ee.Image(10.0).pow(img.select('..').divide(10.0)))
}
function toDB(img) {
  return img.select([]).addBands(ee.Image(img).log10().multiply(10.0));
}

// Implementation by Andreas Vollrath (ESA), inspired by Johannes Reiche (Wageningen)
function terrainCorrection(image) {
  var imgGeom = image.geometry();
  var srtm = ee.Image("USGS/3DEP/10m").clip(imgGeom); // 10m dem
  // to natural
  var sigma0Pow = ee.Image.constant(10).pow(image.divide(10.0));

  // Article ( numbers relate to chapters)
  // 2.1.1 Radar geometry
  var theta_i = image.select('angle');
  var phi_i = ee.Terrain.aspect(theta_i)
    .reduceRegion(ee.Reducer.mean(), theta_i.get('system:footprint'), 1000)
    .get('aspect');

  // 2.1.2 Terrain geometry
  var alpha_s = ee.Terrain.slope(srtm).select('slope');
  var phi_s = ee.Terrain.aspect(srtm).select('aspect');

  // 2.1.3 Model geometry
  // reduce to 3 angle
  var phi_r = ee.Image.constant(phi_i).subtract(phi_s);

  // convert all to radians
  var phi_rRad = phi_r.multiply(Math.PI / 180);
  var alpha_sRad = alpha_s.multiply(Math.PI / 180);
  var theta_iRad = theta_i.multiply(Math.PI / 180);
  var ninetyRad = ee.Image.constant(90).multiply(Math.PI / 180);

  // slope steepness in range (eq. 2)
  var alpha_r = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan();

  // slope steepness in azimuth (eq 3)
  var alpha_az = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan();

  // local incidence angle (eq. 4)
  var theta_lia = (alpha_az.cos().multiply((theta_iRad.subtract(alpha_r)).cos())).acos();
  var theta_liaDeg = theta_lia.multiply(180 / Math.PI);
  // 2.2
  // Gamma_nought_flat
  var gamma0 = sigma0Pow.divide(theta_iRad.cos());
  var gamma0dB = ee.Image.constant(10).multiply(gamma0.log10());
  var ratio_1 = gamma0dB.select('VV').subtract(gamma0dB.select('VH'));

  // Volumetric Model
  var nominator = (ninetyRad.subtract(theta_iRad).add(alpha_r)).tan();
  var denominator = (ninetyRad.subtract(theta_iRad)).tan();
  var volModel = (nominator.divide(denominator)).abs();

  // apply model
  var gamma0_Volume = gamma0.divide(volModel);
  var gamma0_VolumeDB = ee.Image.constant(10).multiply(gamma0_Volume.log10());

  // we add a layover/shadow mask to the original implmentation
  // layover, where slope > radar viewing angle
  var alpha_rDeg = alpha_r.multiply(180 / Math.PI);
  var layover = alpha_rDeg.lt(theta_i);

  // shadow where LIA > 90
  var shadow = theta_liaDeg.lt(85);

  // // calculate the ratio for RGB vis
  // var ratio = gamma0_VolumeDB.select('VV').subtract(gamma0_VolumeDB.select('VH'));

  // var output = gamma0_VolumeDB.addBands(ratio).addBands(alpha_r).addBands(phi_s).addBands(theta_iRad)
  //   .addBands(layover).addBands(shadow).addBands(gamma0dB).addBands(ratio_1);
  var output = image.select('angle').addBands(gamma0_VolumeDB.select(['VV', 'VH'])).updateMask(layover).updateMask(shadow)

  return output;//image.addBands(output.select(['VV', 'VH'], ['VV', 'VH']),null,true);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Refined Lee Speckle - by Guido Lemoine
// Coded as in SNAP 3.0 S1TBX https://goo.gl/cJmzvj
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function RefinedLee(img) {
  // img must be in natural units, i.e. not in dB!
  // Set up 3x3 kernels
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

  var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel);
  var sample_var = variance3.neighborhoodToBands(sample_kernel);

  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());

  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);

  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);

  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));

  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);

  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());

  //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
  //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  // print('DEBUG sigmaV', sigmaV);

  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0],
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

  // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

  // and add the bands for rotated kernels
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }

  // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());

  // A finally generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

  var b = varX.divide(dir_var);

  var result = dir_mean.add(b.multiply(img.subtract(dir_mean)));

  // print('DEBUG result', result);

  // Convert the array image to a scaler image.
  result = result.arrayGet(0);
  // print('DEBUG result (after conversion)', result);

  return img.select([]).addBands(result);
}


//============= Function: logistic function fitting for SG filter ================
function SG_filter(col, time_interval,window,order){
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
  // var ddiff = ee.Date(img.get('date')).difference(ee.Date(start_date), 'hour')
  var ddiff = ee.Number(img.get(time_interval))
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
  var ddiff = ee.Number(img.get(time_interval))
  return img.addBands(ee.Image(1).toFloat().rename('constant'))
      .addBands(ee.Image(ddiff).toFloat().rename('t'))
      .addBands(ee.Image(ddiff).pow(ee.Image(2)).toFloat().rename('t2'))
      .updateMask(img.mask().reduce('min'))
  }).sort('date')
  }


  // For the remainder, use s1res as a list of images
  // set time stamp for regression and sorting

  var col_list = col_res.toList(col_res.size())
  var runLength = ee.List.sequence(0, col_res.size().subtract(window))//.aside(print,'runLength')
  // Run the SG solver over the series, and return the smoothed image version
  var sg_series = runLength.map(function(i) {
    var ref = ee.Image(col_list.get(ee.Number(i).add(half_window)))
    var h = ref.get('date')
    // Convert the collection to an array.
    var img_list = bands.map(function(b){

      var array_ = col_res.select(ee.List([b]).cat(indepSelectors)).toArray();


     // Get a slice corresponding to the window_size of the SG smoother
      var subarray = array_.arraySlice(imageAxis, ee.Number(i).int(), ee.Number(i).add(window).int())

      var array_mask = col_res.select(ee.List([b])).map(function(i){ return i.mask()}).toArray();
      var subarray_mask = array_mask.arraySlice(imageAxis, ee.Number(i).int(), ee.Number(i).add(window).int())
      .arraySlice(bandAxis, 0, 1)
      // make sure there is more validated pixels in windows (subarray) than variables(orders)
      var mask_overdetermined = subarray_mask.arrayReduce('sum', [imageAxis])
      .gt(ee.ImageCollection([ee.Image(order+3)]).toArray())
      .arrayRepeat(bandAxis,  nbands + order + 1)
      subarray = subarray.arrayMask(mask_overdetermined)//.unmask(ee.ImageCollection([ee.Image(0)]).toArray())

      var predictors = subarray.arraySlice(bandAxis, nbands, nbands + order + 1);
      var response = subarray.arraySlice(bandAxis, 0, nbands); //
      var coeff = predictors.matrixSolve(response)

      coeff = coeff.arrayProject([0]).arrayFlatten(coeffFlattener)

      var sged = coeff.multiply(ref.select(indepSelectors)).reduce(ee.Reducer.sum()).rename([ee.String(b).cat('_sg')])
      return ref.select([b]).addBands(sged)
    })
    if(bands.size().gt(1)){
    var img_sged = ee.Image(img_list.slice(1,bands.size()).iterate(function(i,f){ return ee.Image(f).addBands(i)},img_list.get(0)))
    }else{ var img_sged = ee.Image(img_list.get(0))}
    return ref.select([]).addBands(img_sged)
  })

  return ee.ImageCollection(sg_series)
}


// per-band gamma filter
function GammaMap(image) {
  // Cf. https://github.com/senbox-org/s1tbx/blob/master/s1tbx-op-sar-processing/src/main/java/org/esa/s1tbx/sar/gpf/filtering/SpeckleFilters/GammaMap.java
  // which implements Lopes et al, IGARSS 1990, 2409-2412.
  // See: https://www.researchgate.net/publication/224270891_Maximum_A_Posteriori_Speckle_Filtering_And_First_Order_Texture_Models_In_Sar_Images.
  // This is the equivalent of the getGammaMapValue() method
  var enl = 4
  var ksize = 5

  // Convert image from dB to natural values
  var nat_img = ee.Image(10.0).pow(image.divide(10.0));

  // Square kernel, ksize should be odd (typically 3, 5 or 7)
  var weights = ee.List.repeat(ee.List.repeat(1,ksize),ksize);

  // ~~(ksize/2) does integer division in JavaScript
  var kernel = ee.Kernel.fixed(ksize,ksize, weights, ~~(ksize/2), ~~(ksize/2), false);

  // Get mean and variance
  var mean = nat_img.reduceNeighborhood(ee.Reducer.mean(), kernel);
  var variance = nat_img.reduceNeighborhood(ee.Reducer.variance(), kernel);

  // "Pure speckle" threshold
  var ci = variance.sqrt().divide(mean);  // square root of inverse of enl

  // If ci <= cu, the kernel lies in a "pure speckle" area -> return simple mean
  var cu = 1.0/Math.sqrt(enl)

  // If cu < ci < cmax the kernel lies in the low textured speckle area -> return the filtered value
  var cmax = Math.sqrt(2.0) * cu

  var alpha = ee.Image(1.0 + cu*cu).divide(ci.multiply(ci).subtract(cu*cu));
  var b = alpha.subtract(enl + 1.0)
  var d = mean.multiply(mean).multiply(b).multiply(b).add(alpha.multiply(mean).multiply(nat_img).multiply(4.0*enl));
  var f = b.multiply(mean).add(d.sqrt()).divide(alpha.multiply(2.0));

  // If ci > cmax do not filter at all (i.e. we don't do anything, other then masking)

  // Compose a 3 band image with the mean filtered "pure speckle", the "low textured" filtered and the unfiltered portions
  var bands = toDB(mean.updateMask(ci.lte(cu)).addBands(toDB(f.updateMask(ci.gt(cu)).updateMask(ci.lt(cmax)))).addBands(image.updateMask(ci.gte(cmax))))
  var composed_img = ee.Image(bands).select('VV.*').reduce(ee.Reducer.sum()).rename('VV').addBands(ee.Image(bands).select('VH.*').reduce(ee.Reducer.sum()).rename('VH'))
  return image.select([]).addBands(composed_img);
}
