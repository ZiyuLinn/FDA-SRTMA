from initial import *

# %%
# sensitivity test

band_lst = inputs.bandNames().getInfo()
isRS = True
version = 'MEIpca'
type_lst = ['DBF', 'ENF','DNF','nonF','W']
type_num = len(type_lst)
band_type_lst = ['SMA','STMA','RTMA','SRTMA']
SITE_NAME_model_lst = ['UNDE','CHEQ','STEI'] #

SITE_NAME_lst = SITE_NAME_model_lst
scale_lst = [10]
# CLUSTER_LST = [3,5,7,9,11]
CLUSTER_LST = [1]
R=20

for SITE_NAME_m in SITE_NAME_model_lst:
    for sample_type in range(3):

        library_lst= []

        for SITE_NAME_m_ in SITE_NAME_m.split('-'):

            pureEMs_name = glob.glob(r'Y:\Unmixing_Western_GreatLake\V2\\purepixels_{}_{}_R{}_thres*.csv'.format(SITE_NAME_m,R))[0]
            library_2d = pd.read_csv(pureEMs_name, index_col=0)
            if sample_type != 'all':
                library_2d = library_2d[library_2d.tile_num %3 != sample_type]


            library_2d['type'] = library_2d[['DNF', 'ENF', 'DBF', 'nonF']].apply(
                lambda d: [c for c, v in zip(['DNF', 'ENF', 'DBF', 'nonF'], d.values) if v == d.max()][0], axis=1)
            library_2d = library_2d[library_2d[['DNF', 'ENF', 'DBF', 'nonF']].max(axis=1)==100]

            VV_bands = [c for c in library_2d.columns if c.__contains__('VV_')]
            B8_bands = [c for c in library_2d.columns if c.__contains__('B8_')]

            library_2d['type'] = library_2d.apply(
                lambda d: 'W' if (np.mean(d[B8_bands].values) < 580) else  d['type'],
                axis=1)

            groundtruth = ee.Image(
                'users/linziyu/WesternGreatLake/groundtruth_{}_clip_10m'.format(SITE_NAME_m_))
            bound_m = groundtruth.geometry().buffer(2000).bounds()


            if type_lst.__contains__('W'):
                url = inputs.select(band_lst).addBands(water_mask) \
                    .getDownloadURL(
                    {'scale': 100, 'crs': crs, 'format': 'NPY', 'region': bound_m, 'name': "water endmembers"})
                print(url)

                response = requests.get(url)
                data = np.load(io.BytesIO(response.content))
                arr = np.array([list(i) for i in data.flatten() if i[0] > 0])
                arr = arr.reshape(arr.shape[0], len(arr[0]))
                df_w = pd.DataFrame(data=arr,
                                    columns=band_lst+['type'])
                df_w = df_w[df_w.type == 1]
                df_w['type'] = 'W'
                library_2d =library_2d.append(df_w)


        SITE_NAME = SITE_NAME_m

        # nonlocal EM
        # for SITE_NAME in [s for s in ['CHEQ','UNDE','STEI'] if s != SITE_NAME_m]:
        for scale in scale_lst:
            groundtruth = ee.Image(
                'users/linziyu/WesternGreatLake/groundtruth_{}_clip_10m'.format(SITE_NAME))
            bound = groundtruth.geometry().bounds()
            sub_tile = bound.coveringGrid('EPSG:4326', 5000)
            tileNum = sub_tile.size().getInfo()
            crs = groundtruth.projection()

            input_img = inputs.clip(bound).setDefaultProjection(crs, None, 10).reduceResolution(ee.Reducer.mean(), False, scale * scale / 100 + 1) \
                .reproject(crs, None, scale)

            input_img = input_img.addBands([ee.Image(input_img.select('.*' + b + '_.*')).reduce('mean').rename(b) for b in bands_dict['SMA']])

            for band_type in band_type_lst:
                # bands = np.unique([k for j in bands_dict[band_type] for k in band_lst if k.__contains__(j) ]).tolist()
                # ================create vector list ========================
                vector_lst = []
                xbar_lst = []

                feature_lst = [band_type]
                bands_f_lst = []
                for feature_type_ in feature_lst:
                    if band_type == 'SMA':
                        bands_ = bands_dict[feature_type_]
                        for b in bands_:
                            c_selected = [c for c in library_2d.columns if c.__contains__(b+'_')]
                            library_2d[b] = library_2d[c_selected].mean(axis=1)
                    else:
                        bands_ = np.unique([k  for j in bands_dict[feature_type_] for k in band_lst if k.__contains__(j)]).tolist()
                    data_denoise_lst = []

                    for i in range(type_num):
                        type_ = type_lst[i]
                        type_arr = library_2d[library_2d['type'] == type_][bands_].values
                        data_denoise_lst.append(type_arr)
                    if feature_type_[0] == 'R':
                        fdim_num = type_num - 2
                    else:
                        fdim_num = type_num - 1
                    clf = LDA(solver='svd')
                    clf.fit(library_2d[bands_].values, library_2d.type.map({"DBF": 0, "DNF": 1, 'ENF': 2, 'nonF': 3, 'W': 4
                                                          }).values)
                    xbar_ = clf.priors_ @ clf.means_
                    vector_ = clf.scalings_[:,0:fdim_num]

                    vector_lst.append(vector_)
                    xbar_lst.append(xbar_)
                    bands_f_ = [feature_type_+str(f) for f in range(fdim_num)]
                    bands_f_lst.append(bands_f_)
                    # X_new = clf.transform(library_2d[bands_].values).astype(np.float64)
                    X_new = (library_2d[bands_].values-xbar_)@vector_
                    df_fish = pd.DataFrame(data=X_new,
                                           columns=bands_f_)

                    library_2d[bands_f_] =   df_fish[bands_f_].values
                bands_f = np.concatenate(bands_f_lst)


                for n_cluster in CLUSTER_LST:
                    print('============cluster {}============'.format(n_cluster))

                    # ============== endmember selection method =================
                    EM_lst = []
                    EM_f_lst = []
                    EM_num_lst = []
                    EM_num_f_lst = []
                    data_denoise_lst = []

                    for type_ in type_lst:
                        type_arr = library_2d[library_2d['type'] == type_][bands_].values
                        type_f_arr = library_2d[library_2d['type'] == type_][bands_f].values


                        if type_arr.shape[0] > n_cluster:
                            kmeans_ = KMeans(n_clusters=n_cluster, random_state=2022).fit(type_arr)
                            kmeans_f = KMeans(n_clusters=n_cluster, random_state=2022).fit(type_f_arr)
                            EM_lst.append(kmeans_.cluster_centers_)
                            EM_f_lst.append(kmeans_f.cluster_centers_)
                            EM_num_lst.append(n_cluster)
                            EM_num_f_lst.append(n_cluster)
                        else:
                            EM_lst.append(type_arr)
                            EM_f_lst.append(type_f_arr)
                            EM_num_lst.append(type_arr.shape[0])
                            EM_num_f_lst.append(type_arr.shape[0])


                    EM_lst = [list(l) for l in np.concatenate(EM_lst)]
                    EM_f_lst = [list(l) for l in np.concatenate(EM_f_lst)]


                    tile_lst = ee.List([])
                    for tile_id in range(tileNum):
                        print("tile_id:{}".format(tile_id))

                        tile_bound = ee.Feature(sub_tile.toList(tileNum).get(tile_id)).geometry().buffer(10).bounds()
                        tile_img = input_img.select(bands_).clip(tile_bound)#ee.Kernel.square(window_size)
                        tile_mask = tile_img.mask().reduce('min')
                        tile_img = tile_img.updateMask(tile_mask)
                        result_tile = tile_img.select([]).set('features', band_type).set(
                                'cluster_num', n_cluster).set('resolution', scale)



                        print('-------------------{}--{}--{}m--em:{}--{}--cluster{}-------------'.format(SITE_NAME,
                                                                           SITE_NAME_m,scale,sample_type,band_type, n_cluster))
                        fisher_img = ee.Image(tile_img).select([])
                        for n,feature_type_ in enumerate(feature_lst):

                            if feature_type_[0] == 'R':
                                fdim_num = type_num - 2
                            else:
                                fdim_num = type_num - 1
                            if band_type == 'SMA':
                                bands_ = bands_dict[feature_type_]
                            else:
                                bands_ = np.unique([k for j in bands_dict[feature_type_] for k in band_lst if k.__contains__(j)]).tolist()
                            vector_ = vector_lst[n]
                            xbar_ = xbar_lst[n]

                            xbar_img = ee.Image(xbar_.T.tolist())
                            arrayImg2D = ee.Image(tile_img.select(bands_).subtract(xbar_img)).toArray().toArray(1)
                            vector_img = ee.Image(ee.Array(vector_.T.tolist()))
                            # fisher_img_ = vector_img.matrixMultiply(arrayImg2D) \
                            fisher_img_ = vector_img.matrixMultiply(arrayImg2D) \
                                .arrayProject([0]).arrayFlatten([['dim_' +feature_type_+'_'+ str(d) for d in range(fdim_num)]])
                            fisher_img = fisher_img.addBands(fisher_img_)

                        fisher_band_lst = fisher_img.bandNames()#.getInfo()


                        # # =============== GEE unmix() function =================
                        if n_cluster<=15:
                            if n_cluster==1:
                                num = 3
                            else:
                                num = 2 #for saving computational time
                            # unmixed_img_MESMA =MESMA_GEE_v2(EM_lst, EM_num_f_lst, tile_img.select(bands_).toFloat(), bands_, type_lst, EM_num=num)
                            unmixed_img = FCLS_GEE(EM_lst, EM_num_lst, tile_img.select(bands_).toFloat(), bands_, type_lst)
                            unmixed_f_img_MESMA = MESMA_GEE_v2(EM_f_lst, EM_num_f_lst, fisher_img, fisher_band_lst,
                                                          type_lst, EM_num=num)
                            # result_tile = tile_img.select([]).addBands([unmixed_img_MESMA,unmixed_f_img_MESMA])
                            # result_tile = unmixed_img_MESMA.addBands([unmixed_f_img_MESMA]).clipToBoundsAndScale(geometry=tile_bound,scale=scale)
                            result_tile = unmixed_img.addBands([unmixed_f_img_MESMA]).clipToBoundsAndScale(geometry=tile_bound,scale=scale)
                        else:
                            unmixed_img = FCLS_GEE(EM_lst, EM_num_lst, tile_img.select(bands_).toFloat(), bands_, type_lst)
                            unmixed_f_img = FCLS_GEE(EM_f_lst, EM_num_lst, fisher_img, fisher_band_lst,type_lst)
                            result_tile = tile_img.select([]).addBands([unmixed_img,unmixed_f_img])

                        # result_tile = tile_img.select([]).addBands([unmixed_img,unmixed_img_MESMA,unmixed_f_img,unmixed_f_img_MESMA])
                        result_tile = result_tile.updateMask(result_tile.mask().eq(1)).clip(tile_bound) # make sure the edge (fractions) are masked
                        # tile_lst = tile_lst.add(result_tile)
                        ee.batch.Export.image.toAsset(image=result_tile,
                                                      description='MEI_s{}_{}_em{}_R{}_S{}_tile{}'.format(SITE_NAME, band_type,
                                                                                                  n_cluster, R,
                                                                                                  sample_type,tile_id),
                                                      assetId='projects/ee-project1-unmixing/assets/unmixed_results_Revision1/MEI_s{}_{}_em{}_R{}_S{}_tile{}'
                                                      .format(SITE_NAME, band_type, n_cluster, R, sample_type,tile_id),
                                                      region=bound.bounds(), scale=scale, crs=crs, maxPixels=1e13,
                                                      shardSize=64).start()

                    # export image should not be greater than 10485760 bytes (10.48576 mb)
                    # do not use mosaic(); use "mean" will be faster(less computation on vectors/topographic issues)
                    # result_mosaic = ee.ImageCollection(tile_lst).mean().toUint16().clipToBoundsAndScale(geometry=bound.bounds(),scale=scale)
                    #     .updateMask( ee.ImageCollection(tile_lst.slice(50,tile_lst.size())).mean())\
                    #
                    # ee.batch.Export.image.toAsset(image=result_mosaic,
                    #                               description='MEI_s{}_{}_em{}_R{}_{}'.format(SITE_NAME,band_type, n_cluster,R,sample_type),
                    #                               assetId='projects/ee-project1-unmixing/assets/unmixed_results_Revision1/MEI_s{}_{}_em{}_R{}_{}'
                    #                               .format(SITE_NAME,band_type, n_cluster,R,sample_type),
                    #                               region=bound.bounds(), scale=scale, crs=crs,  maxPixels=1e13).start()







