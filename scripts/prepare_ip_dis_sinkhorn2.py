
def data_prepare(count_path,gene_name_list,min_spot=5):
    ### filter genes
    count_df = pd.read_csv(count_path, sep = '\t',index_col = 0)
    select_index = [gene in gene_name_list for gene in list(count_df.index)]
    count_df = count_df.loc[select_index,:]
    select_index = list((count_df != 0).sum(axis=1) >= min_spot)
    count_df = count_df.loc[select_index,:]

    return count_df


def generate_gene_distribution(gene, count_df):
    '''
    input:
        - gene : gene (type : str)(e.g. 'gene', 'gene1+gene2')
        - count_df : normalized count matrix (gene*spot) (type : pd.DataFrame)
    output:
        - flag : {0,1} : 0 -> count_df 中有该gene；1 --》 count_df 无该gene
        - gene_df : Non-zero gene expression value dataframe with spot name as its index (e.g. 'spot1')
                (type : pd.DataFrame)
                
    spot num < min_spot --> trimed
    '''
    
    spot_list = list(count_df.columns)
    
    # check single/multi genes
    if '+' not in gene:
        # simple
        if gene in list(count_df.index):
            gene_df = count_df.loc[gene,:]
            gene_df = gene_df.loc[gene_df > 0]
            if gene_df.shape[0] == 0:
                return 1, None
            return 0, gene_df
        else:
#             print('{} not in counts df'.format(gene))
            return 1, None
    else:
        # complex
        gene = gene.split('+')
        for g in gene:
            if g not in list(count_df.index):
#                 print('{} not in counts df'.format(g))
                return 1, None
        gene_df = count_df.loc[gene[0],:]
        sub_num = len(gene)
        for g in gene[1:]:
            gene_df = gene_df*count_df.loc[g,:]
        gene_df = gene_df.apply(lambda x: np.power(x, 1/sub_num))
        gene_df = gene_df.loc[gene_df > 0]
        if gene_df.shape[0] == 0:
            return 1, None
        return 0, gene_df



def cal_gene_distance(ip, count_df, spot_pos_df, reg = 0.001, iternum = 200, shuffle = False):
    '''
    input:
        - gene_1/2 : genes to be calculated
        - spot_dis_df : spot distance matrix, returned by generate_distance_matrix(spot_loc_path)
        - celltype_percent_1/2 : celltype percentage matrix calculated by spatialDWLS for celltype 
    output:
        - gene_distance : 如果有基因不在count里面返回的都是None
        
    '''
    gene_a, gene_b = ip.split(' - ')
    gene_a = gene_a.strip().replace('(','').replace(')','')
    gene_b = gene_b.strip().replace('(','').replace(')','')
    
    flag_a, gene_df_a = generate_gene_distribution(gene_a,count_df)
    flag_b, gene_df_b = generate_gene_distribution(gene_b,count_df)
    
    if flag_a + flag_b > 0:
        return None

    if shuffle == False:
        gene_pos_a = np.array(spot_pos_df.loc[[spot in list(gene_df_a.index) for spot in list(spot_pos_df.index)],['row','col']])
        gene_pos_b = np.array(spot_pos_df.loc[[spot in list(gene_df_b.index) for spot in list(spot_pos_df.index)],['row','col']])
    else:
        joined_spot_list = list(spot_pos_df.index)
        
        r_select_spot_a = random.sample(joined_spot_list, gene_df_a.shape[0])
        r_select_spot_b = random.sample(joined_spot_list, gene_df_b.shape[0])
        gene_pos_a = np.array(spot_pos_df.loc[[spot in r_select_spot_a for spot in list(spot_pos_df.index)],['row','col']])
        gene_pos_b = np.array(spot_pos_df.loc[[spot in r_select_spot_b for spot in list(spot_pos_df.index)],['row','col']])

    cost_matrix = ot.dist(gene_pos_a,gene_pos_b,metric='euclidean')
    gene_distance_dir = ot.sinkhorn2(np.array(gene_df_a),np.array(gene_df_b), cost_matrix/cost_matrix.max(), reg = reg, numItermax=iternum)[0]
    
    
    ### the R-L dir
    cost_matrix = ot.dist(gene_pos_b,gene_pos_a,metric='euclidean')
    gene_distance_rev = ot.sinkhorn2(np.array(gene_df_b),np.array(gene_df_a), cost_matrix/cost_matrix.max(), reg = reg, numItermax=iternum)[0]
    
    gene_distance = (gene_distance_dir+gene_distance_rev)/2

    return gene_distance


def compare_gene_distance(ip,count_df,spot_pos_df,reg = 0.001, iternum = 200, shuffle_num = 500):
    '''
    比较原先的和置换后的distance
    
    main function
    '''
    dis_ori = cal_gene_distance(ip,count_df,spot_pos_df,reg=reg, iternum=iternum)
    
    if dis_ori == None:
#         print('{} None'.format(ip))
        return [ip, dis_ori, None, None, None]
    
    ## shuffle distance 先用avg去比较吧
    dis_shu = []
    for n in range(shuffle_num):
        tmp_d = cal_gene_distance(ip, count_df, spot_pos_df, reg=reg, iternum=iternum, shuffle=True)
        if tmp_d != None:
            dis_shu.append(tmp_d)
    dis_shu_avg = np.average(dis_shu)
    
    dis_r = dis_ori/dis_shu_avg
    
    n_count = 0
    if dis_r > 1:
        for d in dis_shu:
            if dis_ori < d:
                n_count += 1
    else:
        for d in dis_shu:
            if dis_ori > d:
                n_count += 1
    p_val = n_count/shuffle_num
    
#     print('{} ok'.format(ip))
    
    return [ip, dis_ori, dis_shu_avg, dis_r, p_val]


def processing(idx_path,
               count_path,
               output_path,
               cc_ip_df,
               cc_gene_name,
               min_spot=5,
               reg=0.001,
                    iternum=500,
                    shuffle_num=500
              ):

    ## spot postion & spatial counts
    cc_ip_dic = {}
    for i in range(cc_ip_df.shape[0]):
        cc_ip_dic.update({cc_ip_df.iloc[i,2]:{'pathway':cc_ip_df.iloc[i,0],'ip_type':cc_ip_df.iloc[i,1]}})

    pos_df = pd.read_csv(idx_path,sep = '\t',index_col = 0)
    count_df = data_prepare(count_path,cc_gene_name, min_spot=min_spot)

    
    # filter avi ips (all l/r have expression in st data)
    def select_expressed_lr(df_line, all_genes):
        '''
        find out whether both ligand and recepter genes can be found in sc data
        '''
        flag_l = df_line['ligand'] in all_genes
        if df_line['receptor'].count(',') == 0: # receptor with no subunit
            flag_r = df_line['receptor'] in all_genes
        else:  # receptor with multi-subunits
            flag_r = True
            for tmp_r in df_line['receptor'].split(','):
                tmp_flag = tmp_r in all_genes
                flag_r = flag_r and tmp_flag
        return flag_l and flag_r

    all_genes = list(count_df.index)
    select_index = cc_ip_df.apply(lambda x: select_expressed_lr(x, all_genes), axis = 1)
    cc_ip_df = cc_ip_df.loc[select_index,:]
    cc_ip_list = list(cc_ip_df['interaction_name_2'])

    
    pool = ThreadPoolExecutor(4)

    compare_gene_distance_partial = partial(compare_gene_distance, count_df = count_df, spot_pos_df = pos_df, reg=reg, iternum=iternum, shuffle_num=shuffle_num)

    results_iter = pool.map(compare_gene_distance_partial, cc_ip_list)
    res_list_all = [res for res in results_iter if res.count(None) == 0]
    res_df = pd.DataFrame(res_list_all)
    res_df.columns = ['ip','d_ori','d_shu','d_rat','p_val']
    res_df = res_df.sort_values(by='p_val')
    res_df['pathway'] = res_df.apply(lambda x : cc_ip_dic[x[0]]['pathway'], axis=1)
    res_df['ip_type'] = res_df.apply(lambda x : cc_ip_dic[x[0]]['ip_type'], axis=1)
    
    res_df.to_csv('{}/ip_distance_all.tsv'.format(output_path),sep = '\t', index = None)
    
    print('{} done'.format(idx_path.split('/')[-1]))
    

    
if __name__ == '__main__':
    
    import numpy as np
    import ot
    import pandas as pd
    import matplotlib
    from matplotlib import pyplot as plt
    import random
    from functools import partial
    from concurrent.futures import ThreadPoolExecutor
    import os
    import getopt
    import sys

    import warnings
    warnings.filterwarnings('ignore')

    
    optlist,alist=getopt.getopt(sys.argv[1:],'hc:p:o:')
    for opt in optlist:
        if opt[0] == '-h':sys.exit(usage)
        elif opt[0] == '-c':st_count_path = opt[1]
        elif opt[0] == '-p':st_coord_path = opt[1]
        elif opt[0] == '-o':output_path = opt[1]

    
    
    cc_ip_path = '/fs/home/liuzhaoyang/data/cc_ip/cc_ip_all_multi_split_deduplicates.tsv'

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    cc_ip_df = pd.read_csv(cc_ip_path,sep='\t',index_col=0)
    
    all_ips = list(cc_ip_df.index)
    cc_gene_name = []
    for ip in all_ips:
        cc_gene_name += ip.split('_')
    cc_gene_name = list(set(cc_gene_name))

    processing(st_coord_path,st_count_path,output_path,cc_ip_df,cc_gene_name,shuffle_num=500, reg=0.001)

