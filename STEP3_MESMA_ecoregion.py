

from initial import *

# %%
band_lst = inputs.bandNames().getInfo()
print(inputs.bandNames().getInfo())

version = 'EM-LDA3'
type_lst = ['DBF', 'ENF','DNF','nonF','W']
type_num = len(type_lst)
band_type_lst = ['RTMA','STMA','SRTMA']
SITE_NAME_model_lst = ['CHEQ-UNDE-STEI']

scale_lst = [10]
CLUSTER_LST = [1]

for SITE_NAME_m in SITE_NAME_model_lst:

        library_lst = []
        for SITE_NAME_m_ in SITE_NAME_m.split('-'):
            print(SITE_NAME_m_)
            library_2d = \
                pd.read_csv(r'Y:\Unmixing_Western_GreatLake\V2\{}_purepixels_phenoall_SI800_30m.csv'.format(SITE_NAME_m_
                                                                                                    ),
                            index_col=0)


            library_2d['type'] = library_2d[['DNF', 'ENF', 'DBF', 'nonF']].apply(
                lambda d: [c for c, v in zip(['DNF', 'ENF', 'DBF', 'nonF'], d.values) if v == d.max()][0], axis=1)
            library_2d = library_2d[library_2d[['DNF', 'ENF', 'DBF', 'nonF']].max(axis=1)>90]
            library_2d = library_2d[library_2d.SI<800]

            VV_bands = [c for c in library_2d.columns if c.__contains__('VV_')]
            B8_bands = [c for c in library_2d.columns if c.__contains__('B8_')]

            library_2d['type'] = library_2d.apply(
                lambda d: 'W' if (np.mean(d[B8_bands].values) < 580) else  d['type'],
                axis=1)


            library_lst.append(library_2d)

        library_2d = pd.concat(library_lst)
        library_2d = library_2d[band_lst + ['type']].dropna()


        for band_type in band_type_lst:
            bands = np.unique([k for j in bands_dict[band_type] for k in band_lst if k.__contains__(j) ]).tolist()
            # ================create vector list ==========================================================
            vector_lst = []
            xbar_lst = []
            feature_lst = [band_type]
            bands_f_lst = []
            for feature_type_ in feature_lst:
                bands_ = np.unique([k  for j in bands_dict[feature_type_] for k in bands if k.__contains__(j)]).tolist()
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
                    type_arr = library_2d[library_2d['type'] == type_][bands].values
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
                for tile_id in range(n_tiles):
                    tile_bound = ee.Feature(bound_list.get(tile_id)).geometry().buffer(50).bounds()
                    tile_img = inputs.select(bands).clip(tile_bound)#ee.Kernel.square(window_size)
                    tile_mask = tile_img.mask().reduce('min')
                    tile_img = tile_img.updateMask(tile_mask)
                    result_tile = tile_img.select([]).set('features', band_type).set(
                            'cluster_num', n_cluster).set('resolution', scale)
                    print(tile_img.bandNames().getInfo())


                    print('---------------------{}m---{}--cluster{}-------------'.format(scale,band_type,n_cluster))
                    fisher_img = ee.Image(tile_img).select([])
                    for n,feature_type_ in enumerate(feature_lst):
                        if feature_type_[0] == 'R':
                            fdim_num = type_num - 2
                        else:
                            fdim_num = type_num - 1
                        bands_ = np.unique([k for j in bands_dict[feature_type_] for k in bands if k.__contains__(j)]).tolist()
                        vector_ = vector_lst[n]
                        xbar_ = xbar_lst[n]

                        xbar_img = ee.Image(xbar_.T.tolist())
                        arrayImg2D = ee.Image(tile_img.select(bands_).subtract(xbar_img)).toArray().toArray(1)
                        vector_img = ee.Image(ee.Array(vector_.T.tolist()))
                        # fisher_img_ = vector_img.matrixMultiply(arrayImg2D) \
                        fisher_img_ = vector_img.matrixMultiply(arrayImg2D) \
                            .arrayProject([0]).arrayFlatten([['dim_' +feature_type_+'_'+ str(d) for d in range(fdim_num)]])
                        fisher_img = fisher_img.addBands(fisher_img_)
                        print(fisher_img.bandNames().getInfo())


                    fisher_band_lst = fisher_img.bandNames().getInfo()


                    # # =============== GEE unmix() function =================
                    if n_cluster<=5:
                        if n_cluster==1:
                            num = 3
                        else:
                            num = 2
                        unmixed_img_MESMA =MESMA_GEE_v2(EM_lst, EM_num_f_lst, tile_img.select(bands).toFloat(), bands, type_lst, EM_num=num)
                        unmixed_f_img_MESMA = MESMA_GEE_v2(EM_f_lst, EM_num_f_lst, fisher_img, fisher_band_lst,
                                                      type_lst, EM_num=num)
                        result_tile = tile_img.select([]).addBands([unmixed_img_MESMA,unmixed_f_img_MESMA])
                    else:
                        unmixed_img = FCLS_GEE(EM_lst, EM_num_lst, tile_img.select(bands).toFloat(), bands, type_lst)
                        unmixed_f_img = FCLS_GEE(EM_f_lst, EM_num_lst, fisher_img, fisher_band_lst,type_lst)
                        result_tile = tile_img.select([]).addBands([unmixed_img,unmixed_f_img])

                    # result_tile = tile_img.select([]).addBands([unmixed_img,unmixed_img_MESMA,unmixed_f_img,unmixed_f_img_MESMA])
                    result_tile = result_tile.updateMask(result_tile.mask().eq(1)) # make sure the edge (fractions) are masked
                    tile_lst = tile_lst.add(result_tile)

                    ee.batch.Export.image.toAsset(ee.ImageCollection(tile_lst).mosaic(),
                        'Ecoregion_{}_{}m_em{}_{}_{}'.format(band_type, scale, n_cluster,version,tile_id),
                        'projects/ee-project1-unmixing/assets/unmixed_ecoregion/Ecoregion_{}_{}m_em{}_{}_{}'.format(
                                                                                            band_type,scale,n_cluster,version,tile_id),
                        None, None, tile_bound, scale, crs, None, 1e13).start()







