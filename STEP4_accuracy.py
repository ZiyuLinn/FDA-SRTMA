# -*- encoding: utf-8 -*-
'''
@Contact :   linziyu1996@gmail.com
@License :   (C)Copyright 2021-2026
 
@Modify Time      @Author    @Version   
------------      -------    --------    
4/27/2022 10:03 PM   Ziyu LIN      1.0         

 @Desciption
 ------------------------------------
'''
from initial import *


usgsLandForm = ee.Image("CSP/ERGo/1_0/US/landforms")
usgsDEM = ee.Image("USGS/3DEP/10m")
S2 = ee.ImageCollection("COPERNICUS/S2_SR")

lst = ['DBF', 'DNF', 'ENF', 'nonF']
crs = 'EPSE:4326'


def map_gt(g, scale=10):
    groundtruth = g.expression('(b(0)<=2)?0:b(0)')
    # print(groundtruth.projection().getInfo())
    crs = groundtruth.projection().crs()
    groundtruth = groundtruth.unmask(0).toFloat()
    DBF_count = groundtruth.eq(4).toFloat()
    DNF_count = groundtruth.eq(5).toFloat()
    ENF_count = groundtruth.eq(3).toFloat()
    OV_count = groundtruth.lte(2).toFloat()
    img_ground = DBF_count.addBands([DNF_count, ENF_count, OV_count]) \
        .setDefaultProjection(crs, None, 1).reproject(crs, None, 1) \
        .reduceResolution(ee.Reducer.mean(), False, scale * scale + 1).reproject(crs, None, scale) \
        .rename(lst)

    print(groundtruth.projection().getInfo())
    return img_ground


def print_acc(img, img_ground, scale=20, straitifying=False):
    acc_lst = []
    for i in range(4):
        # print(lst[i])
        acc = accuracy(img.select('.*' + lst[i] + '.*'), img_ground.select(lst[i]), scale, straitifying)
        print('{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{}'.format(acc[0],acc[1],acc[2],acc[3],acc[4]))
        acc_lst.append(acc)

    return acc_lst


def resample_bicubic(img):
    return ee.Image(img).resample('bicubic')


def accuracy(img1_, img2_, scale=20, straitifying=False):
    crs = 'EPSE:4326'
    mask = img1_.mask().Or(img2_.mask())
    acc_scale = scale
    crs = img1_.projection().crs()
    img2 = img2_.updateMask(
        mask)  # .setDefaultProjection(crs, None, 10).reduceResolution(ee.Reducer.mean(), False, scale * scale/100+1).reproject(crs, None,scale)  # .aside(Map.addLayer, {}, 'clip')
    img1 = img1_.updateMask(
        mask)  # .setDefaultProjection(crs, None, 10).reduceResolution(ee.Reducer.mean(), False, scale * scale/100+1).reproject(crs, None,scale)
    # print('1',img1.projection().getInfo())
    # print('2',img2.projection().getInfo())
    bound = img2.geometry()
    img1_b = img1.bandNames().get(0)
    img2_b = img2.bandNames().get(0)
    if straitifying == False:
        img1_mean = ee.Number(
            img1.reduceRegion(ee.Reducer.mean(), bound, acc_scale, None, None, False, 1e9, 16).get(img1_b))
        img2_mean = ee.Number(
            img2.reduceRegion(ee.Reducer.mean(), bound, acc_scale, None, None, False, 1e9, 16).get(img2_b))
        r2_RSS = ee.Number(
            img1.subtract(img2).pow(2).reduceRegion(ee.Reducer.sum(), bound, acc_scale, None, None, False, 1e9, 16).get(
                img1_b))
        r2_TSS = ee.Number(
            img1.subtract(img1_mean).pow(2).reduceRegion(ee.Reducer.sum(), bound, acc_scale, None, None, False, 1e9,
                                                         16).get(img1_b))
        # input_imgicient_of_determination
        R2 = ee.Number(1).subtract(r2_RSS.divide(r2_TSS))  # .aside(print)
        corr_1 = ee.Number(img1.subtract(img1_mean).multiply(img2.subtract(img2_mean))
                           .reduceRegion(ee.Reducer.sum(), bound, acc_scale, None, None, False, 1e9, 16).get(
            img1_b)).abs()
        corr_2 = ee.Number(img1.subtract(img1_mean).pow(2)
                           .reduceRegion(ee.Reducer.sum(), bound, acc_scale, None, None, False, 1e9, 16).get(img1_b))
        corr_3 = ee.Number(img2.subtract(img2_mean).pow(2)
                           .reduceRegion(ee.Reducer.sum(), bound, acc_scale, None, None, False, 1e9, 16).get(img2_b))
        # correlation_input_imgicient
        corr = corr_1.divide(corr_2.sqrt().multiply(corr_3.sqrt()))  # .aside(print)
        # RMSE
        mask = img1.gte(0).updateMask(img2.gte(0))
        RMSE = ee.Number(
            img1.subtract(img2).pow(2).updateMask(mask).reduceRegion(ee.Reducer.mean(), bound, acc_scale, None, None,
                                                                     False,
                                                                     1e9, 16).get(img1_b)).sqrt()  # .aside(print)
        MAE = ee.Number(
            img1.subtract(img2).divide(img2).abs().updateMask(mask).reduceRegion(ee.Reducer.mean(), bound, acc_scale,
                                                                                 None, None,
                                                                                 False,
                                                                                 1e9, 16).get(img1_b))  # .aside(print)
        SLOPE = ee.Number(
            img1.addBands(img2).updateMask(mask).reduceRegion(ee.Reducer.linearFit(), bound, acc_scale,
                                                                                 None, None,False,
                                                                                 1e9, 16).get('scale'))
        weight = -999
    elif straitifying == True:
        straitify_img = img2.divide(0.051).toInt().rename('straitifying')
        acc_lst = []
        sample_lst = []

        # for random in [2,0,23]:
        for random in [2,0,23]:
            sampling = img1.addBands([img2, straitify_img]) \
                .stratifiedSample(100, 'straitifying', bound, scale, img1.projection(), random)
            # print(sampling.getInfo())
            # sampling_img1 = np.array(sampling.aggregate_array(img1_b).getInfo())
            # sampling_img2 = np.array(sampling.aggregate_array(img2_b).getInfo())
            sampling_img = np.array(ee.List([sampling.aggregate_array(img1_b),sampling.aggregate_array(img2_b),sampling.aggregate_array('straitifying')]).getInfo())
            sample_lst.append(sampling_img)
            # print(sampling_img.shape)
            # sample1_lst.append(sampling_img1)
            # sample2_lst.append(sampling_img2)

        sampling_img = np.concatenate(sample_lst,axis=1)
        # print(sampling_img.shape)
        sample1_lst = []
        sample2_lst = []
        for layer in range(0,20):
            sampling_img1 = sampling_img[0,sampling_img[2,:]==layer][0:100]
            sampling_img2 = sampling_img[1,sampling_img[2,:]==layer][0:100]
            # print(layer,len(sampling_img1),len(sampling_img2))
            sample1_lst.append(sampling_img1)
            sample2_lst.append(sampling_img2)
        sampling_img1 = np.concatenate(sample1_lst)
        sampling_img2 = np.concatenate(sample2_lst)
        # sampling_img1 = sampling_img1[0:1000]
        # sampling_img2 = sampling_img2[0:1000]
        R2 = metrics.r2_score(sampling_img1, sampling_img2)
        corr = np.corrcoef(sampling_img1, sampling_img2)[0][1]
        # If squared=True returns MSE value, if False returns RMSE value.
        RMSE = metrics.mean_squared_error(sampling_img1, sampling_img2, squared=False)
        MAE = metrics.mean_absolute_error(sampling_img1, sampling_img2)
        SLOPE  = np.polyfit(sampling_img1,sampling_img2,1)[0]
        weight = len(sampling_img2)
        # print(R2,corr,RMSE)
    acc_1 = [ee.Number(R2).getInfo(), ee.Number(RMSE).getInfo(), ee.Number(corr).pow(2).getInfo(),ee.Number(SLOPE).getInfo(),
            ee.Number(MAE).getInfo(), ee.Number(weight).getInfo()]

    # return np.mean(acc_lst,axis=0)
    return acc_1

# %% accuracy -- kmeans Sensitivity
scale_ = 90
isSf = False
crs = 'EPSG:4326'

em_lst = [1]
scale_lst = [10]

feature_lst = ['SRTMA']

region_lst = ['double']
imageCollection = ee.ImageCollection("projects/ee-project1-unmixing/assets/Kmeans_sensitivity2")  # k-means num

df_lst = []
for site in ['UNDE', 'CHEQ', 'STEI']:
    model_lst = [site]
    img_ground = ee.Image("users/linziyu/WesternGreatLake/groundtruth_{}_clip_10m".format(site))
    crs = img_ground.projection()
    bound = img_ground.geometry().bounds()


    for n_cluster in range(1, 13):
        subtile = bound
        band_type = 'SRTMA'
        version = 'MEI_s{}_{}_em{}_R20'.format(site, band_type, n_cluster)

        for scale_ in [90]:

            img_ground_ = img_ground.setDefaultProjection(crs, None, 10).reproject(crs, None, 10) \
                .reduceResolution(ee.Reducer.mean(), False, scale_ * scale_ / 100 + 1).reproject(crs, None,
                                                                                                 scale_).clip(
                subtile)

            if imageCollection.filterMetadata('system:index', 'contains', version).size().getInfo() == 0:
                continue
            result = imageCollection.filterMetadata('system:index', 'contains', version).first().clip(subtile)

            print(result.get('system:index').getInfo())

            result = result.select('[^W^n]*').addBands(
                result.select(['.*_W', '.*_nonF']).reduce('sum').rename('new_nonF')) \
                .addBands(
                result.select(['.*_W_1', '.*_nonF_1']).reduce('sum').rename('new_nonF_1'))\
                .addBands(
                    result.select(['.*_W_1_1', '.*_nonF_1_1']).reduce('sum').rename('new_nonF_1_1'))

            result = result.divide(10000).setDefaultProjection(crs, None, 10) \
                .reduceResolution(ee.Reducer.mean(), False, scale_ * scale_ / 100 + 1).reproject(crs, None,
                                                                                                 scale_)  # .aside(Map.addLayer, {}, 'clip')

            print('------MESMA--------')
            acc = print_acc(result.select(['.*_DBF', '.*_DNF', '.*_ENF', '.*_nonF'])
                            # .mask(result.select(model+'_rmse_f').unmask(0).lt(100))
                            , img_ground=img_ground_, scale=scale_, straitifying=isSf)
            df = pd.DataFrame(columns=['R2', 'RMSE', 'r2', 'Slope', 'MAE', 'count'], data=acc)
            df['scale'] = scale_
            df['name'] = [n + version for n in ['DBF', 'DNF', 'ENF', 'nonF']]
            print('------fisher MESMA--------')
            acc_f = print_acc(result.select(['.*_DBF_1', '.*_DNF_1', '.*_ENF_1', '.*_nonF_1'])
                              # .mask(result.select(model+'_rmse_f').unmask(0).lt(100))
                              , img_ground=img_ground_, scale=scale_, straitifying=isSf)
            df_f = pd.DataFrame(columns=['R2', 'RMSE', 'r2', 'Slope', 'MAE', 'count'], data=acc_f)
            df_f['scale'] = scale_
            df_f['name'] = [n + version + '_FDA2' for n in ['DBF', 'DNF', 'ENF', 'nonF']]

            acc_f = print_acc(result.select(['.*_DBF_1_1', '.*_DNF_1_1', '.*_ENF_1_1', '.*_nonF_1_1'])
                              # .mask(result.select(model+'_rmse_f').unmask(0).lt(100))
                              , img_ground=img_ground_, scale=scale_, straitifying=isSf)
            df_f2 = pd.DataFrame(columns=['R2', 'RMSE', 'r2', 'Slope', 'MAE', 'count'], data=acc_f)
            df_f2['scale'] = scale_
            df_f2['name'] = [n + version + '_FDA' for n in ['DBF', 'DNF', 'ENF', 'nonF']]
            df_lst.append(df)
            df_lst.append(df_f)
            df_lst.append(df_f2)

pd.concat(df_lst).to_csv(r'Y:\Unmixing_Western_GreatLake\accuracy_{}m_{}.csv'.format(scale_,'Kmeans_sensitivity4'))



# %% TODO accuracy-witin site
scale_ = 90
R = 20
isSf = True
crs = 'EPSG:4326'

em_lst = [1]
scale_lst = [10]
EM_num = 1
feature_lst =['SRTMA']
# feature_lst = ['RTMA','STMA','SRTMA']


imageCollection = ee.ImageCollection("projects/ee-project1-unmixing/assets/unmixed_results_Revision1") #k-means num

slope = ee.Terrain.slope(ee.Image("USGS/3DEP/10m"))
dem = ee.Image("USGS/3DEP/10m")
df_lst = []
for site in ['UNDE','CHEQ','STEI']:

    model_lst = [site]
    version_lst = ['MEI_s{}_{}_em{}_R{}_S{}'.format(site,f,EM_num,R, n)  for f in feature_lst for n in range(3)]

    img_ground = ee.Image("users/linziyu/WesternGreatLake/groundtruth_{}_clip_10m".format(site))
    crs = img_ground.projection()
    bound = img_ground.geometry().bounds()
    boundtiles = bound.coveringGrid('EPSG:4326', 2000)
    size = boundtiles.size().getInfo()

    for k, version in enumerate(version_lst):
        print(version)
        region = version.split('_')[-1]
        if region=='all':
            subtile = bound
        elif region=='double':
            region = 0 if region=='single' else 1
            subtile = ee.FeatureCollection(
                ee.List([r for r in range(size) if r % 2 != region]).map(lambda n: boundtiles.toList(size + 1).get(n)))
        else:
            region = int(region[1:])
            subtile = ee.FeatureCollection(
                ee.List([r for r in range(size) if r%3==region]).map(lambda n: boundtiles.toList(size + 1).get(n)))

        # for scale_ in [30,50,70,90,110]:
        for scale_ in [90]:
            img_ground_ = img_ground.setDefaultProjection(crs, None, 10).reproject(crs, None, 10) \
                .reduceResolution(ee.Reducer.mean(), False, scale_ * scale_ / 100 + 1).reproject(crs, None, scale_).clip(
                subtile)

            if imageCollection.filterMetadata('system:index', 'contains', version).size().getInfo()==0:
                continue
            result = imageCollection.filterMetadata('system:index', 'contains', version).mosaic().clip(subtile)

            print(result.get('system:index').getInfo())
            # continue


            result = result.select('[^W^n]*').addBands(
                result.select(['.*_W', '.*_nonF']).reduce('sum').rename('new_nonF')) \
                .addBands(
                result.select(['.*_W_1', '.*_nonF_1']).reduce('sum').rename('new_nonF_1'))


            result = result.divide(10000).setDefaultProjection(crs, None, 10) \
                .reduceResolution(ee.Reducer.mean(), False,scale_*scale_/100+1).reproject(crs, None,
                                                                      scale_)  # .aside(Map.addLayer, {}, 'clip')

            print('------MESMA--------')
            acc = print_acc(result.select(['.*_DBF', '.*_DNF', '.*_ENF', '.*_nonF'])
                               , img_ground=img_ground_, scale=scale_, straitifying=isSf)
            df = pd.DataFrame(columns=['R2', 'RMSE', 'r2', 'Slope','MAE', 'count'], data=acc)
            df['scale']= scale_
            df['name'] = [n +  version for n in ['DBF', 'DNF', 'ENF', 'nonF']]
            print('------fisher MESMA--------')
            acc_f = print_acc(result.select(['.*_DBF_1', '.*_DNF_1', '.*_ENF_1', '.*_nonF_1'])
                               # .mask(result.select(model+'_rmse_f').unmask(0).lt(100))
                               , img_ground=img_ground_, scale=scale_, straitifying=isSf)
            df_f = pd.DataFrame(columns=['R2', 'RMSE', 'r2', 'Slope', 'MAE','count'], data=acc_f)
            df_f['scale']= scale_
            df_f['name'] = [n+version+'_FDA' for n in ['DBF', 'DNF', 'ENF', 'nonF']]
            df_lst.append(df)
            df_lst.append(df_f)
pd.concat(df_lst).to_csv(r'Y:\Unmixing_Western_GreatLake\Accuracy_MEIpca_{}m_R{}_EM{}_{}_SMA.csv'.format(scale_,R,EM_num,'within-site'))

  # %% TODO accuracy-cross site
scale_ = 90
isSf = True
crs = 'EPSG:4326'
EM_num = 1
R=20
feature_lst = ['RTMA', 'STMA', 'SRTMA']

imageCollection = ee.ImageCollection("projects/ee-project1-unmixing/assets/unmixed_results_Revision1")  # k-means num
slope = ee.Terrain.slope(ee.Image("USGS/3DEP/10m"))
dem = ee.Image("USGS/3DEP/10m")
df_lst = []
for site in ['UNDE','CHEQ','STEI']:
    model_lst = [s for s in ['UNDE','STEI','CHEQ'] if s !=site]
    version_lst = ['MEI_s{}_m{}-SAM_{}'.format(site,m,f) for f in feature_lst for m in model_lst]

    img_ground = ee.Image("users/linziyu/WesternGreatLake/groundtruth_{}_clip_10m".format(site))
    crs = img_ground.projection()
    bound = img_ground.geometry().bounds()
    boundtiles = bound.coveringGrid('EPSG:4326', 2000)
    size = boundtiles.size().getInfo()

    for k, version in enumerate(version_lst):
        print(version)
        region = version.split('_')[-1][1:]
        region = 'all'
        if region == 'all':
            subtile = bound
        elif region == 'double':
            region = 0 if region == 'single' else 1
            subtile = ee.FeatureCollection(
                ee.List([r for r in range(size) if r % 2 != region]).map(lambda n: boundtiles.toList(size + 1).get(n)))
        else:
            region = int(region)
            subtile = ee.FeatureCollection(
                ee.List([r for r in range(size) if r % 3 != region]).map(lambda n: boundtiles.toList(size + 1).get(n)))

        # for scale_ in [30,50,70,90,110]:
        for scale_ in [scale_]:
            img_ground_ = img_ground.setDefaultProjection(crs, None, 10).reproject(crs, None, 10) \
                .reduceResolution(ee.Reducer.mean(), False, scale_ * scale_ / 100 + 1).reproject(crs, None,
                                                                                                 scale_).clip(
                subtile)

            if imageCollection.filterMetadata('system:index', 'contains', version).size().getInfo() == 0:
                continue
            result = imageCollection.filterMetadata('system:index', 'contains', version).mean().clip(subtile)

            print(result.get('system:index').getInfo())
            # continue

            result = result.select('[^W^n]*').addBands(
                result.select(['.*_W', '.*_nonF']).reduce('sum').rename('new_nonF')) \
                .addBands(
                result.select(['.*_W_1', '.*_nonF_1']).reduce('sum').rename('new_nonF_1'))

            result = result.divide(10000).setDefaultProjection(crs, None, 10) \
                .reduceResolution(ee.Reducer.mean(), False, scale_ * scale_ / 100 + 1).reproject(crs, None,
                                                                                                 scale_)  # .aside(Map.addLayer, {}, 'clip')

            # print('------MESMA--------')
            acc = print_acc(result.select(['.*_DBF', '.*_DNF', '.*_ENF', '.*_nonF'])
                            # .mask(result.select(model+'_rmse_f').unmask(0).lt(100))
                            , img_ground=img_ground_, scale=scale_, straitifying=isSf)
            df = pd.DataFrame(columns=['R2', 'RMSE', 'r2', 'Slope', 'MAE', 'count'], data=acc)
            df['scale'] = scale_
            df['name'] = [n + version for n in ['DBF', 'DNF', 'ENF', 'nonF']]
            print('------fisher MESMA--------')
            acc_f = print_acc(result.select(['.*_DBF_1', '.*_DNF_1', '.*_ENF_1', '.*_nonF_1'])
                              # .mask(result.select(model+'_rmse_f').unmask(0).lt(100))
                              , img_ground=img_ground_, scale=scale_, straitifying=isSf)
            df_f = pd.DataFrame(columns=['R2', 'RMSE', 'r2', 'Slope', 'MAE', 'count'], data=acc_f)
            df_f['scale'] = scale_
            df_f['name'] = [n + version + '_FDA' for n in ['DBF', 'DNF', 'ENF', 'nonF']]
            df_lst.append(df)
            df_lst.append(df_f)
pd.concat(df_lst).to_csv(r'Y:\Unmixing_Western_GreatLake\Accuracy_MEIpca_{}m_R{}_EM{}_{}_SAM_mean.csv'.format(scale_, R,EM_num,'cross-site'))
