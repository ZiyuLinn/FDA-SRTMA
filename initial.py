# -*- encoding: utf-8 -*-
'''
@File    :   F-Unmixing_createEndMembers.py
@Contact :   linziyu1996@gmail.com
@License :   (C)Copyright 2021-2026

@Modify Time      @Author    @Version
------------      -------    --------
4/25/2022 12:01 PM   Ziyu LIN      1.0

 ------------------------------------

'''

#%%
# Trigger the authentication flow for Earth Engine
# !earthengine --ee_config EE_CONFIG authenticate

# init the ee object
import ee
ee.Initialize()

# import all needed package
import numpy as np
import pandas as pd
import time,os,io,sys,glob,requests,itertools
import matplotlib.pyplot as plt
# import mesma
# from mesma.core.mesma import MesmaModels,MesmaCore
import sklearn.metrics as metrics
import sklearn.metrics as metrics
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.cluster import KMeans,OPTICS,DBSCAN, AgglomerativeClustering
# %matplotlib inline
import matplotlib.pyplot as plt
#%%=======================global parometers===========================
projection = 'EPSG:4326'
crs = 'EPSG:4326'
scale = 10
pi = 3.141593

ecoregion = ee.FeatureCollection("users/linziyu/WesternGreatLake/northern_highland_level4")
boundtiles = ecoregion.geometry().bounds().coveringGrid('EPSG:4326', 20000).filterBounds(ecoregion)
n_tiles = boundtiles.size().getInfo()
bound_list = boundtiles.toList(n_tiles)

s2_img = ee.ImageCollection('projects/ee-project1-unmixing/assets/S2_ecoregion')\
    .select('month.*sg').mosaic()
s1_img = ee.ImageCollection('projects/ee-project1-unmixing/assets/S1_ecoregion')
# fill nan and filter outliners
s1_img_mask = s1_img.select('QA').mosaic().gt(30)
s1_img =s1_img.select('month.*sg').map(lambda img: ee.Algorithms.If(ee.Image(img).bandNames().contains('month11_VV'),ee.Image(img),
ee.Image(img).addBands(ee.Image(img).select('.*month4.*').add(ee.Image(img).select('.*month10.*')).divide(2).toInt16()\
.rename(ee.Image(img).select('.*month4.*').bandNames().map(lambda b: ee.String(b).replace('month4','month11'))))))\
.select('month.*sg').mosaic().updateMask(s1_img_mask)

slope = ee.Terrain.slope(ee.Image("USGS/3DEP/10m"))
slope_mask = slope.select('slope').lt(10)
water_mask = s1_img.select(['.*VV_.*']).reduce('mean').lt(445).And(s2_img.select(['.*B8_.*']).reduce('mean').lt(580))
lowbiomass_mask = s1_img.select(['.*VV_.*']).reduce('mean').expression('(b(0)>445)&(b(0)<929)?1:0').updateMask(water_mask.Not())
urban_mask = s1_img.select(['.*VV_.*']).reduce('mean').expression('(b(0)>4366)?1:0')

RF_bands = [ 'B2', 'B3', 'B4','B5','B6','B7', 'B8', 'B8A','B11','B12']
SAR_bands = ['VV','VH']
bands_dict = {'SRTMA': RF_bands+SAR_bands,'RTMA': SAR_bands,'STMA': RF_bands,'SMA': RF_bands}
inputs = s1_img.addBands(s2_img).select('month[1-9].*_sg')
band_lst = inputs.bandNames().getInfo()

# print(inputs.bandNames().getInfo())

# %% Functions

def threshold_otsu(counts, bin_centers, *args, **kwargs) -> float:
    """Find the threshold value for a bimodal histogram using the Otsu method.

    If you have a distribution that is bimodal (AKA with two peaks, with a valley
    between them), then you can use this to find the location of that valley, that
    splits the distribution into two.

    From the SciKit Image threshold_otsu implementation:
    https://github.com/scikit-image/scikit-image/blob/70fa904eee9ef370c824427798302551df57afa1/skimage/filters/thresholding.py#L312
    """
    # counts, bin_edges = np.histogram(x, *args, **kwargs)
    # bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2

    # class probabilities for all possible thresholds
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(counts * bin_centers) / weight1
    mean2 = (np.cumsum((counts * bin_centers)[::-1]) / weight2[::-1])[::-1]

    # Clip ends to align class 1 and class 2 variables:
    # The last value of ``weight1``/``mean1`` should pair with zero values in
    # ``weight2``/``mean2``, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    idx = np.argmax(variance12)
    threshold = bin_centers[idx]
    return threshold,idx

def resample_bicubic(img):
    return ee.Image(img).resample('bicubic')

#============== mesma for earth engine ====================
# create EM list for different combinations
# 2-EM models :
# inputs: EM_lst, tile_img,bands,type_name_lst
# outputs: best_unmixed_img with min_error_bands(in MAE)
# version2 resolve the computational problem
def MESMA_GEE_v2(EM_lst,EM_num_lst,tile_img,bands,type_name_lst,EM_num=2):
    # create endmembers combination
    type_lst = ['{}-{}'.format(j,t) for i,t in enumerate(type_name_lst) for j in range(EM_num_lst[i])]

    numbers = [i for i in range(len(EM_lst))]
    # n-EMs combinations based on  itertools.combinations
    # index_lst_ = []
    # for k in range(2,EM_num+1):
    #
    #     combinations = list(itertools.combinations(numbers, k))
    #     index_lst_ = index_lst_ + combinations
    index_lst_ = list(itertools.combinations(numbers, EM_num))

    index_lst = []
    for inds in index_lst_:
        type_lst_sub = [type_lst[ind] for ind in inds]
        if len(np.unique([t.split('-')[1] for t in type_lst_sub])) == EM_num:
            index_lst.append(inds)

    print(index_lst,len(index_lst))
    N = 50
    if len(index_lst)>N:
        counter_lst = [[c for c in range(i-N,i)] for i in range(N,len(index_lst),N)]
        if len(index_lst)%N>0:
            counter_lst=counter_lst+[[c for c in range(len(index_lst)-len(index_lst)%N,len(index_lst))]]
    else:
        counter_lst = [[c for c in range(len(index_lst))]]

    results_sum_all = tile_img.select([])
    for c_num,counters in enumerate(counter_lst):
        # print(counters)
        error_bands = tile_img.select([])
        unmixed_bands = tile_img.select([])
        for i in counters:
            inds = index_lst[i]
            EM_lst_sub = [EM_lst[ind] for ind in inds]
            type_lst_sub = [type_lst[ind] for ind in inds]
            if len(np.unique([t.split('-')[1] for t in type_lst_sub]))==1:
                continue
            # else:
            #     counter= counter+1
            # # =============== GEE unmix() function =================
            unmixed_img = tile_img.select(bands).toFloat().unmix(EM_lst_sub,True, True)\
                .rename(['frac_{}_{}'.format(i,t) for t in type_lst_sub])
            # print( unmixed_img.bandNames().getInfo())
            # =========== RMSE for SMA============
            tile_img_pred = ee.Image(ee.Array(EM_lst_sub).transpose()).toFloat()\
                          .matrixMultiply(unmixed_img.toArray().toArray(1))\
                          .arrayProject([0]).arrayFlatten([bands])
            # error = tile_img.select(bands)\
            #     .subtract(tile_img_pred).abs().reduce('mean').multiply(1000).toInt16().rename('MAE_{}'.format(i))
            error = tile_img.select(bands).spectralDistance(tile_img_pred, 'sam') \
                .multiply(10000).abs().toInt16().rename('SAM_{}'.format(i))

            error_bands = error_bands.addBands(error) #.addBands([error.rename('frac_{}_{}'.format(i,t)) for t in type_lst_sub])
            unmixed_bands = unmixed_bands.addBands(unmixed_img)

        error_min = error_bands.reduce('min').rename('SAM_'+str(c_num))
        error_min_mask = error_bands.eq(error_min)
        # unmixed_bands = unmixed_bands.updateMask(MAE_min_mask).unmask(0)
        best_unmixed_masked = tile_img.select([])
        for i in counters:
            inds = index_lst[i]
            type_lst_sub = [type_lst[ind] for ind in inds]
            # if len(np.unique([t.split('-')[1] for t in type_lst_sub]))==1:
            #     continue
            best_unmixed_band = unmixed_bands.select(['frac_{}_{}'.format(i,t) for t in type_lst_sub])\
                .where(error_min_mask.select('SAM_'+str(i)).Not(),0)
            best_unmixed_masked = best_unmixed_masked.addBands(best_unmixed_band)


        results_sum = tile_img.select([])
        # print( np.unique([text.split('_')[-1].split('-')[1] for text in unmixed_bands.bandNames().getInfo()]))
        for t in np.unique([text.split('_')[-1].split('-')[1] for text in unmixed_bands.bandNames().getInfo()]):
            results_sum = results_sum.addBands(best_unmixed_masked.select(['.*'+t]).reduce('sum').rename(t+'_'+str(c_num)))
        #     .where(results_sum.gt(1.1),1)
        results_sum_all = results_sum_all.addBands(results_sum.addBands([error_min]))
    # print(results_sum_all.bandNames().getInfo())
    error_min = results_sum_all.select('SAM.*').reduce('min')
    error_min_mask = results_sum_all.select('SAM.*').eq(error_min)
    results_sum_all_masked =tile_img.select([])
    for i in  range(len(counter_lst)):
        s = '_'+str(i)
        results_sum_all_masked = results_sum_all_masked.addBands(results_sum_all.select(['.*'+s]).where(error_min_mask.select('SAM'+s).Not(),0))
    results = tile_img.select([])
    for t in type_name_lst:
        results = results.addBands(
            results_sum_all_masked.select([t+'.*']).reduce('sum').rename(t))

    results = results.divide(results.reduce('sum')).unmask(0).multiply(10000).addBands(error_min.rename('SAM_min')).clip(tile_img.geometry()).toUint16()
    return results

#============== FCLS for earth engine ====================
# inputs: EM_lst, tile_img,bands,type_name_lst
# outputs: best_unmixed_img with min_error_bands
def FCLS_GEE(EM_lst,EM_num_lst,tile_img,bands,type_name_lst):
    # create endmembers combination
    index_lst = []

    unmixed_img = tile_img.select(bands).toFloat().unmix(EM_lst, True, True)\
    .rename(['{}_{}'.format(t,n) for i,t in enumerate(type_name_lst) for n in range(EM_num_lst[i])])
    # =========== RMSE for SMA============
    tile_img_pred = ee.Image(ee.Array(EM_lst).transpose()) \
        .matrixMultiply(unmixed_img.toArray().toArray(1)) \
        .arrayProject([0]).arrayFlatten([bands])
    error = tile_img.select(bands) \
        .subtract(tile_img_pred).abs().reduce('mean').rename('MAE')

    results_sum = tile_img.select([])
    for t in type_name_lst:
        results = unmixed_img.select([t + '.*']).reduce('sum').rename(SITE_NAME_m + '_' + t)
        results_sum = results_sum.addBands(results)
    results_sum = results_sum.divide(results_sum.reduce('sum')).unmask(0).multiply(10000).addBands(error).clip(tile_img.geometry()).toUint16()
    return results_sum


# Function to perform Fisher's transformation on the data
def fisher_trans(data,type_num):
    num = sum([len(data[i]) for i in range(type_num)])
    m_in = [np.mean(data[i], axis=0) for i in range(type_num)]
    # the global mean
    m = np.mean(m_in, axis=0)
    # m = np.sum([sum(data[i]) for i in range(f)],axis=0)/num
    # print(np.array(m_in).shape, np.array(m).shape)

    S_w = sum([np.dot((data[i] - m_in[i]).T, (data[i] - m_in[i])) for i in range(type_num)]) / num
    S_b = sum([np.outer((m_in[i] - m).T, (m_in[i] - m)) * len(data[i]) for i in range(type_num)]) / num

    A = np.dot(np.linalg.inv(S_w), S_b)
    value, vector = np.linalg.eig(A)
    return vector


# Function to draw a 3D plot of the data
def drawplot(data_, w_):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    class_lst = []
    for i, s, c, n in zip(range(len(data_)), ['bs', 'go', 'rp', '*'], ['r', 'b', 'g', 'k'], type_lst):
        class1 = np.dot(data_[i], w_).astype(np.float64) / 1000
        print(class1.shape)
        # class_lst.append(class1)
        # print(class1.shape)
        # cluster = np.dot(data_c[i], w)
        # plt.plot(class1[:, 0], class1[:, 1], s, label=n,alpha=0.02, c=c)

        ax.scatter(class1[:, 0], class1[:, 1], class1[:, 2], alpha=0.02, c=c, label=n)
    # for i,s,c,n in zip(range(len(data_)),['bs','go','rp','*'],['r', 'b', 'g', 'k'],type_lst):
    #     norm0 = np.linalg.norm(np.vstack(class_lst)[:,0])
    #     norm1 = np.linalg.norm(np.vstack(class_lst)[:,1])
    #     print(norm0,norm1)
    #     plt.plot(class_lst[i][:, 0]/norm0, class_lst[i][:, 1]/norm1,s, label=n,alpha=0.02, c=c)

    plt.xlabel('Fisher Transformation 1', fontsize=15)
    plt.ylabel('Fisher Transformation 2', fontsize=15)
    ax.set_zlabel('Fisher Transformation 3', fontsize=15)

    plt.legend(fontsize=15, frameon=False, ncol=4)
    plt.show()