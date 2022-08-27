def trans_ip_symbol(ip):
    ip = ip.split('_')
    if len(ip) > 2:
        part_b = '({})'.format('+'.join(ip[1:]))
    else:
        part_b = ip[1]
    return '{} - {}'.format(ip[0],part_b)

def trans_ip_symbol_rev(ip):
    ip = ip.split(' - ')
    if ip[1].count('+') > 0:
        part_b = '_'.join(ip[1][1:-1].split('+'))
    else:
        part_b = ip[1]
    return '{}_{}'.format(ip[0],part_b)



if __name__ == '__main__':
    
    import pandas as pd
    import pickle as pkl
    import numpy as np

    data_base_dir = './'

    sc_raw_count_path = '/fs/home/liuzhaoyang/project/cci_evaluation/human_intestinal/data/processed/intestinal_sc_count.tsv'
    sc_raw_meta_path = '/fs/home/liuzhaoyang/project/cci_evaluation/human_intestinal/data/processed/intestinal_sc_meta.tsv'

    sc_round1_count_path = '{}/sc_count_round1.tsv'.format(data_base_dir)
    sc_round1_meta_path = '{}/sc_meta.tsv'.format(data_base_dir)

    st_coord_path = '{}/st_coord.tsv'.format(data_base_dir)

    target_ct_distype_path = '{}/target_ct_distype_sr.pkl'.format(data_base_dir)

    target_ctp_ip_path = '{}/target_ct_pair_ip_dic.pkl'.format(data_base_dir)
    target_ctp_ip_tf_path = '{}/target_ct_pair_ip_tf_tg_dic.pkl'.format(data_base_dir)

    d_rat_path = '{}/ip_dis_sinkhorn2/ip_distance_all.tsv'.format(data_base_dir)

    simu_spot_cell_dic_path = '{}/simu_spot_cell_dic.pkl'.format(data_base_dir)

    ################
    # reading data #
    ################
    sc_count_df_raw = pd.read_csv(sc_raw_count_path, sep='\t', index_col=0)
    sc_meta_df_raw = pd.read_csv(sc_raw_meta_path, sep='\t', index_col=0, header=None)
    sc_meta_df_raw.columns = ['celltype']

    sc_count_df_r1 = pd.read_csv(sc_round1_count_path, sep='\t', index_col=0)
    sc_meta_df_r1 = pd.read_csv(sc_round1_meta_path, sep='\t', index_col=0, header=None)
    sc_meta_df_r1.columns = ['celltype']

    st_coord_df_r1 = pd.read_csv(st_coord_path, sep='\t', index_col=0)

    d_rat_df = pd.read_csv(d_rat_path, sep='\t', index_col=0)

    with open(target_ct_distype_path, 'rb') as f:
        target_ct_distype_sr = pkl.load(f)
    with open(target_ctp_ip_path, 'rb') as f:
        target_ctp_ip_dic = pkl.load(f)
    with open(target_ctp_ip_tf_path, 'rb') as f:
        target_ctp_ip_tf_dic = pkl.load(f)
    with open(simu_spot_cell_dic_path, 'rb') as f:
        simu_spot_cell_dic = pkl.load(f)


    ###################
    # modify sc count #
    ###################

    # subset raw sc count df (keep the same gene & cell with simulated round1)
    sc_count_df_enriched = sc_count_df_raw.loc[list(sc_count_df_r1.index),list(sc_count_df_r1.columns)]
    # sc_meta_df_raw = sc_meta_df_raw.loc[list(sc_meta_df_r1.index),:]

    ### select & keep short/long target ips

    # filter short/long-range ips
    top_per = 0.1
    pval_cutoff = 0.01
    max_index = int(d_rat_df.shape[0]*top_per)

    # filtering short-range interactions
    d_rat_df = d_rat_df.sort_values(by=['d_rat','p_val'])
    tmp_d_rat_df = d_rat_df.loc[d_rat_df.p_val <= pval_cutoff,:]
    max_index = max(max_index, tmp_d_rat_df.shape[0])
    short_ips = [trans_ip_symbol_rev(ip) for ip in list(d_rat_df.iloc[:max_index,:].index)]
    short_ips = set(short_ips)

    # filtering long-range interactions
    d_rat_df = d_rat_df.sort_values(by=['d_rat','p_val'],ascending=[False,True])
    tmp_d_rat_df = d_rat_df.loc[d_rat_df.p_val <= pval_cutoff,:]
    max_index = max(max_index, tmp_d_rat_df.shape[0])
    long_ips = [trans_ip_symbol_rev(ip) for ip in list(d_rat_df.iloc[:max_index,:].index)]
    long_ips = set(long_ips)


    # near/far celltypes
    near_cts = list(target_ct_distype_sr.loc[target_ct_distype_sr == 0].index)
    near_cts_rev = ['{}|{}'.format(ctp.split('|')[1], ctp.split('|')[0]) for ctp in near_cts]
    near_cts += near_cts_rev

    far_cts = list(target_ct_distype_sr.loc[target_ct_distype_sr == 1].index)
    far_cts_rev = ['{}|{}'.format(ctp.split('|')[1], ctp.split('|')[0]) for ctp in far_cts]
    far_cts += far_cts_rev


    # filter keep ips and corresponding tf, tg per celltype pair
    print('> filtering ips')
    gene_barcode_dic = {} # gene & enriched barcodes, values are set
    for tmp_ctp in target_ctp_ip_tf_dic:

        tmp_target_ips = list(target_ctp_ip_tf_dic[tmp_ctp].keys())
        if tmp_ctp in near_cts:
            tmp_ct_distype = 'near'
            tmp_avi_ips = set(tmp_target_ips).intersection(short_ips)
        elif tmp_ctp in far_cts:
            tmp_ct_distype = 'far'
            tmp_avi_ips = set(tmp_target_ips).intersection(long_ips)
        print('{} - {} - {} ips kept'.format(tmp_ctp, tmp_ct_distype, len(tmp_avi_ips)))

        del_ips = [ip for ip in tmp_target_ips if ip not in tmp_avi_ips]
        for ip in del_ips:
            del target_ctp_ip_tf_dic[tmp_ctp][ip]

        # link genes & barcodes
        ct_a, ct_b = tmp_ctp.split('|')
        ct_a_barcodes = list(sc_meta_df_r1.loc[sc_meta_df_r1.celltype == ct_a,:].index)
        ct_b_barcodes = list(sc_meta_df_r1.loc[sc_meta_df_r1.celltype == ct_b,:].index)

        for ip in target_ctp_ip_tf_dic[tmp_ctp]:
            ips = ip.split('_')

            # ligand
            lig = ips[0]
            if lig not in gene_barcode_dic:
                gene_barcode_dic.update({
                    lig:set(ct_a_barcodes)
                })
            else:
                for bar in ct_a_barcodes:
                    gene_barcode_dic[lig].add(bar)

            # receptor
            for rec in ips[1:]:
                if rec not in gene_barcode_dic:
                    gene_barcode_dic.update({
                        rec:set(ct_b_barcodes)
                    })
                else:
                    for bar in ct_b_barcodes:
                        gene_barcode_dic[rec].add(bar)

            if len(target_ctp_ip_tf_dic[tmp_ctp][ip]['tf']) > 0:
                # tf
                for tf in target_ctp_ip_tf_dic[tmp_ctp][ip]['tf']:
                    if tf not in gene_barcode_dic:
                        gene_barcode_dic.update({
                            tf:set(ct_b_barcodes)
                        })
                    else:
                        for bar in ct_b_barcodes:
                            gene_barcode_dic[tf].add(bar)

                # tg
                for tg in target_ctp_ip_tf_dic[tmp_ctp][ip]['tg']:
                    if tg not in gene_barcode_dic:
                        gene_barcode_dic.update({
                            tg:set(ct_b_barcodes)
                        })
                    else:
                        for bar in ct_b_barcodes:
                            gene_barcode_dic[tg].add(bar)


    ## copy avi "enriched" signals to raw sc count
    for tmp_g in gene_barcode_dic:
        tmp_bars = gene_barcode_dic[tmp_g]
        sc_count_df_enriched.loc[tmp_g,tmp_bars] = sc_count_df_r1.loc[tmp_g,tmp_bars]



    ##########################
    # mapping cells to spots #
    ##########################
    # fill the simulated st count data
    all_st_spots = list(st_coord_df_r1.index)
    simulated_st_count_df = pd.DataFrame(
        np.zeros((sc_count_df_enriched.shape[0], len(all_st_spots))),
        index = list(sc_count_df_enriched.index),
        columns = all_st_spots
    )

    for tmp_spot in all_st_spots:
        tmp_exp = round(np.mean(
            sc_count_df_enriched.loc[:,simu_spot_cell_dic[tmp_spot]],
            axis = 1
        ))
        simulated_st_count_df.loc[:,tmp_spot] = tmp_exp
    simulated_st_count_df = simulated_st_count_df.astype(int)






    ##########
    # output #
    ##########

    sc_count_df_enriched.to_csv('{}/sc_count.tsv'.format(data_base_dir), sep='\t')

    with open('{}/target_ctp_ip_tf_tg_avi_dic.pkl'.format(data_base_dir), 'wb') as f:
        pkl.dump(target_ctp_ip_tf_dic, f)

    simulated_st_count_df.to_csv('{}/st_count.tsv'.format(data_base_dir), sep='\t')