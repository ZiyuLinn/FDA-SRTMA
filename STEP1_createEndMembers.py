# -*- encoding: utf-8 -*-
'''
@Contact :   linziyu1996@gmail.com
@License :   (C)Copyright 2021-2026

@Modify Time      @Author    @Version
------------      -------    --------
4/25/2022 12:01 PM   Ziyu LIN      1.0

 @Desciption
 ------------------------------------
'''
from initial import *

# %% sampling end menbers for ALL features in 3 groundtruth sites
# quick downloading in tiles and visulization
for site in [ 'UNDE','STEI', 'CHEQ']:
    print('site:' + site)

    CRS = 'EPSG:4326'

    SITE = site
    SCALE = 20

    groudtruth =ee.Image('users/linziyu/WesternGreatLake/groundtruth_' + site + '_clip_10m')
    CRS = groudtruth.projection()
    ground_per = groudtruth.multiply(100)
    BOUND = groudtruth.geometry().bounds()
    boundtiles = BOUND.bounds().coveringGrid('EPSG:4326', 3000)
    crs = groudtruth.projection()
    pure_ground = ground_per.eq(100).reduce('max')

    inputs = s2_img.addBands([s1_img]).toInt().clip(
        ecoregion.geometry().bounds())
    inputs = inputs.updateMask(inputs.mask().reduce('min'))
    band_lst = inputs.bandNames().getInfo()

    for R in range(50, 51, 5):
        #  the MEI map was generated based on ee.Image.spectralGradient() provided by earth engine
        MEI = ee.Image(
            "projects/ee-project1-unmixing/assets/endmember_extraction/MEI_R{}_{}".format(R,SITE))\
            .divide(pi)
        histogram =MEI.multiply(1000).toInt16()\
        .sample(BOUND, SCALE, CRS, None, 2000, 2023, True, 16)\
        .reduceColumns(ee.Reducer.histogram(200, 5), ['AMEE']).get('histogram')
        counts = np.array(ee.Dictionary(histogram).get('histogram').getInfo())
        bucketMeans = np.array(ee.Dictionary(histogram).get('bucketMeans').getInfo())
        # th,idx = threshold_otsu(counts[bucketMeans>700], bucketMeans[bucketMeans>700])
        th,idx = threshold_otsu(counts, bucketMeans)
        print('threshold MEI = 0.{}'.format(str(th)[0:3]))
        plt.bar( bucketMeans,counts,width=5)
        # plt.plot( bucketMeans,counts,lw=1)
        plt.plot( [th,th],[0,counts[idx]],'r--')
        # plt.plot( [th,th],[0,counts[bucketMeans>700][idx]],'r--')
        plt.ylim([0,max(counts)*1.1])
        plt.xlim([500,1000])
        plt.title("MEI_R{}_{}: 0.{}".format(R,SITE,th))
        plt.show()

        # Create purepixel map, overlap with the groundtruth map and then extarct the pure pixels.
        purity = MEI.multiply(1000)
        purity_mask = purity.gt(th)
        mask = purity_mask.And(pure_ground).And(slope_mask).focalMin(1).setDefaultProjection(CRS, None, 10)\
                .reduceResolution(ee.Reducer.mean(), False, SCALE * SCALE / 100 + 1).eq(1)
        inputs_site = inputs.addBands([ground_per, purity]).updateMask(mask).toInt16()
            # .setDefaultProjection(CRS, None, 10)
            # .reduceResolution(ee.Reducer.mean(), False, SCALE * SCALE / 100 + 1) \
            # .reproject(crs, None, SCALE)
        band_lst = inputs_site.bandNames().getInfo()
        band_lst_purity = purity.bandNames().getInfo()

        num_tile = int(boundtiles.size().getInfo())
        df_all = pd.DataFrame(columns=band_lst)
        for t in range(num_tile):
            print('subtile:{}'.format(t))
            subtile = ee.Feature(boundtiles.toList(num_tile + 1).get(t)).geometry()
            url = inputs_site.clip(subtile) \
                .getDownloadURL(
                {'scale': SCALE, 'crs': CRS, 'format': 'NPY', 'region': subtile, 'name': "thumbnail_subtiles"})
            response = requests.get(url)
            data = np.load(io.BytesIO(response.content))
            arr = np.array([list(i) for i in data.flatten() if all([ii>=0 for ii in i])])
            arr = arr.reshape(arr.shape[0], len(arr[0]))


            df = pd.DataFrame(data=arr, columns=band_lst)
            df['tile_num'] = t
            df_all = pd.concat([df_all,df])

        df_all.to_csv(
            r'Y:\Unmixing_Western_GreatLake\V2\\purepixels_{}_R{}_thres{}.csv'.format(site,R, str(th)[0:3]))

