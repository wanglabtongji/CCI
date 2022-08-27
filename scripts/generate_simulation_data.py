import sys
import os
import getopt


usage='Usage:'+sys.argv[0]
usage+='''<Required>[Options]
    <Required>
    --sc_count: sc count matrix path
    --sc_meta: sc meta path
    --st_coord: ST coord path
    --st_meta: ST meta path
    --target_ct: simulated target cell types
    [Options]
    --cc_ip: cci database path
    --rec_tf: receptor-tf database path
    --tf_target: tf-target gene database path
    --output_dir: simulation output dir
    --min_cell_per_spot: min number of cells mapping to a single spot
    --max_cell_per_spot: max number of cells mapping to a single spot
    --min_cells: remove genes expressed in less than min_cells
    --nip_per_ctp: simulated interaction number in each cell type pair
    --add_fc: permutate & raise fold change for each gene
    --n_target_gene: the number of target genes are selected for each LR pair
'''


###############
# set default #
###############

cc_ip = './cci_database/cc_ip_all_multi_split.tsv'
rec_tf = './data_simulation/extdata/new_ligand_receptor_TFs.txt'
tf_target = './data_simulation/extdata/tf_target.txt'

simulate_output_dir = './'

# set the max and min number of cells mapping to a single spot
min_cell_per_spot = 2
max_cell_per_spot = 5
min_cells = 20
nip_per_ctp = 30
target_ip_cell_per = 0
target_ip_cell_per_flag = False

add_sigma = 0.01 # avoid 0 when calculate log2

# permutate & raise fold change for each gene
add_fc = 1

tf_flag = True # whether simulate TFs & target genes
n_target_gene = 10 # how many target genes are selected for each LR pair


if len(sys.argv) < 2 or not sys.argv[1].startswith('-'):
    sys.exit(usage)

optlist,alist=getopt.getopt(sys.argv[1:],
                            '',
                            [
                                'help=',
                                'sc_count=',
                                'sc_meta=',
                                'st_coord=',
                                'st_meta=',
                                'cc_ip=',
                                'rec_tf=',
                                'tf_target=',
                                'target_ct=',
                                'output_dir=',
                                'min_cell_per_spot=',
                                'max_cell_per_spot',
                                'min_cells=',
                                'nip_per_ctp=',
                                'add_fc=',
                                'n_target_gene='
                            ])
for opt, arg in optlist:
    if opt in ['--help']:
        sys.exit(usage)
    elif opt in ['--sc_count']:
        sc_count_path = arg
    elif opt in ['--sc_meta']:
        sc_meta_path = arg
    elif opt in ['--st_coord']:
        st_pos_path = arg
    elif opt in ['--st_meta']:
        st_meta_path = arg
    elif opt in ['--cc_ip']:
        cc_ip_path = arg
    elif opt in ['--rec_tf']:
        rec_tf_path = arg
    elif opt in ['--tf_target']:
        tf_target_path = arg
    elif opt in ['--target_ct']:
        target_ct_path = arg
    elif opt in ['--output_dir']:
        simulate_output_dir = arg
    elif opt in ['--min_cell_per_spot']:
        min_cell_per_spot = int(arg)
    elif opt in ['--max_cell_per_spot']:
        max_cell_per_spot = int(arg)
    elif opt in ['--min_cells']:
        min_cells = int(arg)
    elif opt in ['--nip_per_ctp']:
        nip_per_ctp = int(arg)
    elif opt in ['--add_fc']:
        nip_per_ctp = float(arg)
    elif opt in ['--n_target_gene']:
        nip_per_ctp = int(arg)


def cal_ct_avg_dis(tmp_meta_df, tmp_pos_df, 
                   ct_name_correct_dic=None, 
                   nspot_cutoff=5,
                  meta_col = 'celltype',
                  pos_x_col = 'row',
                  pos_y_col = 'col'):
    
    avg_dis_sr, ct_name_list = [], []
        
    if ct_name_correct_dic != None:
        tmp_meta_df[meta_col] = tmp_meta_df[meta_col].apply(lambda x: ct_name_correct_dic[x] if x in ct_name_correct_dic else x)
       
    ct_list = list(set(list(tmp_meta_df[meta_col])))  
        
    for n, ct_x in enumerate(ct_list):
        select_index = [tmp_ct == ct_x for tmp_ct in list(tmp_meta_df.loc[:,meta_col])]
        ct_x_barcode = list(tmp_meta_df.loc[select_index,:].index)
        
        if len(ct_x_barcode) < nspot_cutoff:
            continue
            
        select_barcode = [bar in ct_x_barcode for bar in list(tmp_pos_df.index)]
        sub_pos_df_x = tmp_pos_df.loc[select_barcode,:]

        for m in range(n+1, len(ct_list)):
            ct_y = ct_list[m]
            select_index = [tmp_ct == ct_y for tmp_ct in list(tmp_meta_df.loc[:,meta_col])]
            ct_y_barcode = list(tmp_meta_df.loc[select_index,:].index)
            
            if len(ct_y_barcode) < nspot_cutoff:
                continue
            
            select_barcode = [bar in ct_y_barcode for bar in list(tmp_pos_df.index)]
            sub_pos_df_y = tmp_pos_df.loc[select_barcode,:]
            
            ct_name_list.append('{}|{}'.format(ct_x, ct_y))
            
            tmp_dis_matrix = ot.dist(np.array(sub_pos_df_x.loc[:,[pos_x_col,pos_y_col]]), 
                                     np.array(sub_pos_df_y.loc[:,[pos_x_col,pos_y_col]]),
                                     metric='euclidean')

            avg_dis = (np.average(np.amin(tmp_dis_matrix,0)) + np.average(np.amin(tmp_dis_matrix,1)))/2

            avg_dis_sr.append(avg_dis)
    avg_dis_sr = pd.Series(avg_dis_sr, index=ct_name_list)
            
    return avg_dis_sr


def generate_ct_distype(avg_dis_sr):
    ct_dis_array = np.array([[d,1] for d in list(avg_dis_sr)])
    estimator = KMeans(n_clusters=3)
    estimator.fit(ct_dis_array)
    label_pred = estimator.labels_
    label_pred = pd.Series(label_pred, index=list(avg_dis_sr.index))

    avg_label_dis = pd.Series([np.average(avg_dis_sr.loc[label_pred == l]) for l in set(label_pred)])
    avg_label_rank = avg_label_dis.rank()
    min_label = list(avg_label_rank[avg_label_rank == 1].index)[0]
    max_label = list(avg_label_rank[avg_label_rank == 3].index)[0]
    mid_label = [l for l in set(label_pred) if l != min_label and l != max_label][0]
    ct_distype_sr = label_pred.replace({min_label:0,max_label:1,mid_label:2})
    print('near(0): {}; far(1): {}; mid(2): {}'.format(ct_distype_sr.loc[ct_distype_sr==0].shape[0], ct_distype_sr.loc[ct_distype_sr==1].shape[0], ct_distype_sr.loc[ct_distype_sr==2].shape[0]))
    
    return ct_distype_sr


def generate_gene_distribution(gene, count_df):
    '''
    input:
        - gene : gene (type : str)(e.g. 'gene', 'gene1+gene2')
        - count_df : normalized count matrix (gene*spot) (type : pd.DataFrame)
    output:
        - flag : {0,1} : 0 -> count_df 中有该gene；1 --》 count_df 无该gene
        - gene_df : Non-zero gene expression value dataframe with spot name as its index (e.g. 'spot1')
                (type : pd.DataFrame)
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

    

def imputate_target_gene(tmp_gene_sr, tmp_target_barcodes, 
                        imputate_cell_per = 1, 
                        r_ncell = 5,
                        adjust_max = 1.2,
                        adjust_min = 0.6,
                        return_barcodes = False):
    
    # non-zero expression
    tmp_non_zero_sr = tmp_gene_sr.loc[tmp_gene_sr > 0]
    nonzero_barcodes = list(tmp_non_zero_sr.index)
    base_exp = np.min(tmp_non_zero_sr)
    
    # imputation needed cells
    tmp_target_cell_sr = tmp_gene_sr.loc[tmp_target_barcodes]
    zero_barcodes = list(tmp_target_cell_sr.loc[tmp_target_cell_sr < 1].index)
    # if needed, only imputate a part of cells
    imputate_barcodes = random.sample(zero_barcodes, 
                                      k=round(len(zero_barcodes)*imputate_cell_per))
    
    # imputate
    for bar in imputate_barcodes:
        # random select non-zero expression & fill with their average
        if len(nonzero_barcodes) > r_ncell:
            tmp_avg_barcodes = random.sample(nonzero_barcodes, k=r_ncell)
        else:
            tmp_avg_barcodes = nonzero_barcodes
        adjust_cor = random.uniform(adjust_min, adjust_max)
        tmp_exp = round(np.mean(tmp_non_zero_sr.loc[tmp_avg_barcodes])*adjust_cor)
        tmp_exp = np.max([tmp_exp, base_exp])
        tmp_gene_sr.loc[bar] = tmp_exp
    
    if return_barcodes:
        return tmp_gene_sr, imputate_barcodes
    else:
        return tmp_gene_sr
    

if __name__ == '__main__':
    
    import pandas as pd
    import numpy as np
    import ot
    import random
    from sklearn.cluster import KMeans
    import itertools
    import pickle as pkl


    ##########################################################################
        
    ########################
    # simulation celltypes #
    ########################
    
    simulate_target_cts = []
    with open(target_ct_path, 'r') as f:
        for line in f.readlines():
            simulate_target_cts.append(line.strip())
       
    ##########################################################################

    ## reading data
    st_pos_df = pd.read_csv(st_pos_path, sep='\t', index_col = 0)
    st_meta_df = pd.read_csv(st_meta_path, sep = '\t', index_col=0, header = None)
    st_meta_df.columns = ['celltype']
    sc_count_df = pd.read_csv(sc_count_path, sep='\t', index_col=0)
    sc_meta_df = pd.read_csv(sc_meta_path, sep='\t', index_col=0, header=None)
    sc_meta_df.columns = ['celltype']

    all_ip_df = pd.read_csv(cc_ip_path, sep='\t', index_col = 0)
    all_ip_df['ligand'] = all_ip_df.apply(lambda x: x[2].split(' - ')[0].strip(), axis = 1)
    all_ip_df['receptor'] = all_ip_df.apply(lambda x: x[2].split(' - ')[1].replace('(','').replace(')','').replace('+',','), axis = 1)
    all_ip_df = all_ip_df.drop_duplicates()


    # check near/far celltypes
    ct_dis_sr = cal_ct_avg_dis(st_meta_df, st_pos_df,nspot_cutoff=5)
    ct_distype_sr = generate_ct_distype(ct_dis_sr)

    # near ct pair
    print('> all near celltypes\n{}'.format(list(ct_distype_sr.loc[ct_distype_sr == 0].index)))

    # far ct pair
    print('> all far celltypes\n{}'.format(list(ct_distype_sr.loc[ct_distype_sr == 1].index)))


    print('> target celltypes\n{}'.format(', '.join(simulate_target_cts)))

    # plot target celltypes, if needed
    # plot_st_celltype(st_pos_df, st_meta_df, target_ct=simulate_target_cts)

    # subset st meta & position data
    select_index = [ct in simulate_target_cts for ct in list(st_meta_df.loc[:,'celltype'])]
    sub_st_meta_df = st_meta_df.loc[select_index,:]
    sub_st_pos_df = st_pos_df.loc[select_index,:]
    print('> {} spots are selected'.format(sub_st_meta_df.shape[0]))

    # --------------------------------------- #
    # decide how many cell need to be selected from sc data
    # based on spot number of each celltype

    # set the max and min number of cells mapping to a single spot
    min_cell_per_spot = min_cell_per_spot
    max_cell_per_spot = max_cell_per_spot

    # get spot number of target celltypes
    target_ct_nspot_dic = {}
    for tmp_ct in simulate_target_cts:
        target_ct_nspot_dic.update({
            tmp_ct: st_meta_df.loc[st_meta_df.celltype == tmp_ct,:].shape[0]
        })
    print('> spot number of target celltypes\n{}'.format(
        '\n'.join(['{}: {},'.format(tmp_ct, target_ct_nspot_dic[tmp_ct]) for tmp_ct in target_ct_nspot_dic]))
         )

    # set selected cell number for each target celltype
    target_ct_ncell_dic = {}
    for tmp_ct in target_ct_nspot_dic:
        target_ct_ncell_dic.update({
            tmp_ct: round(random.randint(10*min_cell_per_spot, 10*max_cell_per_spot)/10 * target_ct_nspot_dic[tmp_ct])
        })
    print('> randomly selected cell number of target celltypes\n{}'.format(
        '\n'.join(['{}: {},'.format(tmp_ct, target_ct_ncell_dic[tmp_ct]) for tmp_ct in target_ct_ncell_dic]))
         )

    # random select target number of cells from sc data
    r_select_barcodes = []
    for n_ct, tmp_ct in enumerate(target_ct_ncell_dic):
        tmp_barcodes = list(sc_meta_df.loc[sc_meta_df.celltype == tmp_ct,:].index)

        if len(tmp_barcodes) < target_ct_ncell_dic[tmp_ct]:
            r_select_barcodes += tmp_barcodes
        else:
            r_select_barcodes += random.sample(tmp_barcodes, k=target_ct_ncell_dic[tmp_ct])

    select_col = [bar in r_select_barcodes for bar in list(sc_count_df.columns)]
    sub_sc_count_df = sc_count_df.loc[:,select_col]

    select_index = [bar in r_select_barcodes for bar in list(sc_meta_df.index)]
    sub_sc_meta_df = sc_meta_df.loc[select_index,:]
    sub_sc_meta_df = sub_sc_meta_df.loc[list(sub_sc_count_df.columns),:] # actually they should already in the same order, but just in case 
    print('> {} cells are randomly selected'.format(sub_sc_meta_df.shape[0]))


    # only foucs on near/far target celltype pairs
    all_possi_target_ct_pairs = ['|'.join(ctp) for ctp in itertools.permutations(simulate_target_cts,2)]

    near_ct_pairs = list(ct_distype_sr.loc[ct_distype_sr == 0].index)
    reverse_tmp_ct_pairs = ['|'.join(list(reversed(tmp_ctp.split('|')))) for tmp_ctp in near_ct_pairs]
    near_ct_pairs += reverse_tmp_ct_pairs
    far_ct_pairs = list(ct_distype_sr.loc[ct_distype_sr == 1].index)
    reverse_tmp_ct_pairs = ['|'.join(list(reversed(tmp_ctp.split('|')))) for tmp_ctp in far_ct_pairs]
    far_ct_pairs += reverse_tmp_ct_pairs

    near_target_ct_pairs = list(set(near_ct_pairs).intersection(set(all_possi_target_ct_pairs)))
    far_target_ct_pairs = list(set(far_ct_pairs).intersection(set(all_possi_target_ct_pairs)))
    all_target_ct_pairs = near_target_ct_pairs + far_target_ct_pairs

    print('> near celltype pairs\n{}'.format('\n'.join(near_target_ct_pairs)))
    print('> far celltype pairs\n{}'.format('\n'.join(far_target_ct_pairs)))
    # all interactions have been loaded at the beginning


    min_cells = min_cells

    # remove genes expressed in less than 10 cells 
    select_index = list((sub_sc_count_df != 0).sum(axis=1) >= min_cells)
    sub_sc_count_df = sub_sc_count_df.loc[select_index,:]
    print('after filtering, {} genes are left'.format(sub_sc_count_df.shape[0]))



    # filter aviailable lr pairs based on sc data
    all_genes = list(sub_sc_count_df.index)

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

    select_index = all_ip_df.apply(lambda x: select_expressed_lr(x, all_genes), axis = 1)
    all_ip_df = all_ip_df.loc[select_index,:]
    print('after filtering, {} interactions are kept'.format(all_ip_df.shape[0]))




    # random select interactions for each celltype pair
    # the same interactions can be existed in both directions ????
    # maybe that will not influence the results too much ?

    # maybe in far celltype pairs, 
    # contact type interactions should not be selected

    nip_per_ctp = nip_per_ctp
    target_ip_cell_per = target_ip_cell_per
    target_ip_cell_per_flag = target_ip_cell_per_flag

    target_ct_pair_ip_dic = {} # key: celltype_pair; value: "enriched" interactions.store randomly selected ips between celltypes (with direction)
    target_ct_pair_ip_tf_tg_dic = {}
    target_ct_lr_gene_dic = {} # key: celltype; value: "enriched" l/r genes in this celltype. store genes which need to be "enriched" in each celltype
    target_gene_ct_dic = {} # key: l/r gene; value: "enriched" celltypes
    all_available_ct_ip_dic = {} # key: celltype; value: {'l': avialable (ligand) ips, 'r': avialable (receptor) ips }
    all_target_ips = [] # all randomly selected ips across all celltypes
    
    tf_flag = tf_flag # whether simulate TFs & target genes

    all_ips = list(all_ip_df.index)
    all_secreted_ips = list(all_ip_df.loc[all_ip_df.annotation == 'Secreted Signaling',:].index)

    # fill all_available_ct_ip_dic
    for tmp_ct in simulate_target_cts:
        tmp_ct_count_df = sub_sc_count_df.loc[:,sub_sc_meta_df.celltype == tmp_ct]
        ncell_max = tmp_ct_count_df.shape[1]
        if target_ip_cell_per_flag:
            select_index = list((tmp_ct_count_df != 0).sum(axis=1) >= ncell_max*target_ip_cell_per)
        else:
            select_index = list(sub_sc_count_df.index)
        tmp_ct_avi_genes = list(tmp_ct_count_df.loc[select_index,:].index)

        # filter avil ligand/receptor direction interactions
        select_index = [tmp_g in tmp_ct_avi_genes for tmp_g in list(all_ip_df.loc[:,'ligand'])]
        all_available_ct_ip_dic.update({
            tmp_ct: {
                'ligand': list(all_ip_df.loc[select_index,:].index)
            }
        })
        tmp_r_ips = []
        for i, tmp_r in enumerate(list(all_ip_df.loc[:,'receptor'])):
            tmp_r = tmp_r.split(',')
            if sum([t_r in tmp_ct_avi_genes for t_r in tmp_r]) == len(tmp_r):
                tmp_r_ips.append(all_ips[i])
        all_available_ct_ip_dic[tmp_ct].update({
            'receptor': tmp_r_ips
        })

    # randomly select ips
    for tmp_ctp in all_target_ct_pairs:
        tmp_ct_a, tmp_ct_b = tmp_ctp.split('|')

        # randomly select interactions
        tmp_avil_ips = list(set(all_available_ct_ip_dic[tmp_ct_a]['ligand']).intersection(set(all_available_ct_ip_dic[tmp_ct_b]['receptor'])))
        if tmp_ctp in far_target_ct_pairs:
            # for far celltype pairs, only select secreted interactions
            tmp_avil_secreted_ips = list(set(tmp_avil_ips).intersection(set(all_secreted_ips)))
            if len(tmp_avil_secreted_ips) <= nip_per_ctp:
                tmp_target_ips = tmp_avil_secreted_ips
            else:
                tmp_target_ips = random.sample(tmp_avil_secreted_ips, nip_per_ctp)
        else:
            if len(tmp_avil_ips) <= nip_per_ctp:
                tmp_target_ips = tmp_avil_ips
            else:
                tmp_target_ips = random.sample(tmp_avil_ips, nip_per_ctp)

        # save select results
        target_ct_pair_ip_dic.update({
            tmp_ctp: tmp_target_ips
        })
        
        if tf_flag:
            for ip in tmp_target_ips:
                if tmp_ctp not in target_ct_pair_ip_tf_tg_dic:
                    target_ct_pair_ip_tf_tg_dic.update({
                        tmp_ctp:{
                            ip: {'tf':[], 'tg':[]}
                        }
                    })
                else:
                    target_ct_pair_ip_tf_tg_dic[tmp_ctp].update({
                        ip:{'tf':[], 'tg':[]}
                    })

        all_target_ips += tmp_target_ips

        for ct_i, tmp_ct_a in enumerate(tmp_ctp.split('|')):
            # extract ligand/receptor genes
            if ct_i == 0: # tmp_ct_a in the sender (ligand) state
                tmp_genes = [tmp_ip.split('_')[ct_i] for tmp_ip in tmp_target_ips]
            else: # tmp_ct_a in the receiver (receptpr) state
                tmp_genes = []
                for tmp_ip in tmp_target_ips:
                    tmp_genes += tmp_ip.split('_')[ct_i:]

            # fill target_ct_lr_gene_dic
            if tmp_ct_a not in target_ct_lr_gene_dic:
                target_ct_lr_gene_dic.update({
                    tmp_ct_a: tmp_genes
                })
            else:
                target_ct_lr_gene_dic[tmp_ct_a] = list(set(target_ct_lr_gene_dic[tmp_ct_a]+tmp_genes))

            # fill target_gene_ct_dic
            for tmp_g in tmp_genes:
                if tmp_g not in target_gene_ct_dic:
                    target_gene_ct_dic.update({
                        tmp_g: [tmp_ct_a]
                    })
                else:
                    if tmp_ct_a not in target_gene_ct_dic[tmp_g]:
                        target_gene_ct_dic[tmp_g].append(tmp_ct_a)




    if tf_flag:
        # read in ligand-recetor-tf-target_gene data
        # source from cellcall
        rec_tf_df = pd.read_csv(rec_tf_path, sep='\t')
        tf_target_df = pd.read_csv(tf_target_path, sep='\t')
        rec_tf_df = rec_tf_df.iloc[:,5:]
        tf_target_df = tf_target_df.iloc[:,:2]


        # only keep target ips
        all_target_ips_df = pd.DataFrame({
            'ip':all_target_ips,
            'ligand':[ip.split('_')[0] for ip in all_target_ips],
            'receptor':[','.join(ip.split('_')[1:]) for ip in all_target_ips] # complexes are existed in rec_tf_df
        })
        all_target_ips_df.drop_duplicates(inplace=True)

        avi_rec_tf_df = pd.merge(all_target_ips_df, rec_tf_df, left_on='ligand', right_on='Ligand_Symbol')
        avi_rec_tf_df = avi_rec_tf_df.loc[:,['Ligand_Symbol','Receptor_Symbol','TF_Symbol']]
        avi_rec_tf_df = pd.merge(all_target_ips_df, avi_rec_tf_df, left_on='receptor',right_on='Receptor_Symbol')

        # only keep expressed tfs
        avi_gene_df = pd.DataFrame({'expressed':list(sub_sc_count_df.index)})
        avi_rec_tf_df = pd.merge(avi_rec_tf_df, avi_gene_df, left_on='TF_Symbol', right_on='expressed')
        avi_rec_tf_df = avi_rec_tf_df.iloc[:,:-1]
        avi_rec_tf_df.drop_duplicates(inplace=True)

        # only keep expressed target genes
        avi_tf_target_df = pd.merge(tf_target_df, avi_gene_df, left_on='Target_Symbol', right_on='expressed')
        avi_tf_target_df = avi_tf_target_df.iloc[:,:-1]
        avi_tf_target_df.drop_duplicates(inplace=True)

        # merge ligand-recrptor-tf-target 
        avi_rec_tf_target_df = pd.merge(avi_rec_tf_df, avi_tf_target_df, left_on='TF_Symbol', right_on='TF_Symbol')
        avi_rec_tf_target_df = avi_rec_tf_target_df.loc[:,['ip','ligand','receptor','TF_Symbol','Target_Symbol']]
        avi_rec_tf_target_df.columns = ['ip','ligand','receptor','tf','target']



    print('> simulation start')

    # log2 transform original count data
    add_sigma = add_sigma # avoid 0 when calculate log2
    sc_count_df_log2 = np.log2(sub_sc_count_df+add_sigma)

    # permutate & raise fold change for each gene
    add_fc = add_fc
    n_target_gene = n_target_gene # how many target genes are selected for each LR pair


    gene_imputate_barcodes_dic = {}
    all_target_genes = list(target_gene_ct_dic.keys())

    for tmp_gene in all_target_genes:
        tmp_gene_sr = sc_count_df_log2.loc[tmp_gene,:]

        # first permutate barcodes to remove original signals
        random.shuffle(tmp_gene_sr)

        # get "enriched" cell barcodes
        tmp_target_ct = target_gene_ct_dic[tmp_gene]
        tmp_select_index = [ct in tmp_target_ct for ct in list(sub_sc_meta_df.loc[:,'celltype'])]
        tmp_target_barcodes = list(sub_sc_meta_df.loc[tmp_select_index,:].index)

        # imputate target gene in selected cells
        tmp_gene_sr, tmp_imputate_barcodes = imputate_target_gene(tmp_gene_sr, tmp_target_barcodes, 
                                                                  imputate_cell_per=0.6, return_barcodes=True)
        if tf_flag:
            gene_imputate_barcodes_dic.update({
                tmp_gene: tmp_imputate_barcodes
            })

        # add "enriched" signal
        tmp_gene_sr.loc[tmp_target_barcodes] += add_fc

        # cover original signal with "enriched" signal
        sc_count_df_log2.loc[tmp_gene,:] = tmp_gene_sr

    if tf_flag:
        # if needed, TFs & target genes of target LR pair will also be "enriched"
        all_tf_avi_ips = list(avi_rec_tf_target_df.loc[:,'ip'])
        tf_target_barcode_dic = {}

        # select "enriched" tfs & target genes & corresponding "enriched" cells
        done_recs = []
        for tmp_ctp in target_ct_pair_ip_dic:
            tmp_tf_avi_ips = list(set(all_tf_avi_ips).intersection(set(target_ct_pair_ip_dic[tmp_ctp])))
            if len(tmp_tf_avi_ips) == 0:
                continue
            else:
                tmp_recs = [','.join(ip.split('_')[1:]) for ip in tmp_tf_avi_ips]
                for ti, tmp_r in enumerate(tmp_recs):
                    if tmp_r not in done_recs:

                        # get imputation barcodes (receptor gene imputated)
                        if tmp_r.count(',') == 0: # receptor is not a complex
                            tmp_imp_bars = gene_imputate_barcodes_dic[tmp_r]
                        else: # receptor is a complex
                            tmp_imp_bars = []
                            for tr in tmp_r.split(','):
                                tmp_imp_bars += gene_imputate_barcodes_dic[tr]
                            tmp_imp_bars = list(set(tmp_imp_bars))

                        # select TFs & barcodes to be "enriched"
                        tmp_avi_target_df = avi_rec_tf_target_df.loc[avi_rec_tf_target_df.receptor == tmp_r,:]
                        if tmp_avi_target_df.shape[0] > n_target_gene:
                            r_index = random.sample(
                                range(tmp_avi_target_df.shape[0]),
                                k=n_target_gene
                            )
                            tmp_target_genes = list(tmp_avi_target_df.iloc[r_index,-1])
                            tmp_tf_genes = list(set(tmp_avi_target_df.iloc[r_index,-2]))
                        else:
                            tmp_target_genes = list(tmp_avi_target_df.loc[:,'target'])
                            tmp_tf_genes = list(tmp_avi_target_df.loc[:,'tf'])

                        target_ct_pair_ip_tf_tg_dic[tmp_ctp][tmp_tf_avi_ips[ti]]['tf'] = tmp_tf_genes
                        target_ct_pair_ip_tf_tg_dic[tmp_ctp][tmp_tf_avi_ips[ti]]['tg'] = tmp_target_genes
                        
                        # save select tfs & target genes
                        for tg in tmp_target_genes + tmp_tf_genes:
                            if tg not in tf_target_barcode_dic:
                                tf_target_barcode_dic.update({
                                    tg: tmp_imp_bars
                                })
                            else:
                                tf_target_barcode_dic[tg] += tmp_imp_bars
                                tf_target_barcode_dic[tg] = list(set(tf_target_barcode_dic[tg]))

                        # save done receptors
                        done_recs.append(tmp_r)
                        if tmp_r.count(',') > 0:
                            done_recs += tmp_r.split(',')
                    else:
                        continue

        for tmp_gene in tf_target_barcode_dic:
            tmp_tf_targets_sr = sc_count_df_log2.loc[tmp_gene,:]
            random.shuffle(tmp_tf_targets_sr)
            tmp_tf_targets_sr = imputate_target_gene(tmp_tf_targets_sr, tf_target_barcode_dic[tmp_gene])
            tmp_tf_targets_sr.loc[tf_target_barcode_dic[tmp_gene]] += add_fc
            sc_count_df_log2.loc[tmp_gene,:] = tmp_tf_targets_sr


    # transform back
    sc_count_df_enriched = abs(round(np.power(2,sc_count_df_log2)-add_sigma))
    sc_count_df_enriched = sc_count_df_enriched.astype(int)

    print('> simulation done')


    # mapping cells to st data
    all_sc_cells = list(sc_count_df_enriched.columns)
    all_st_spots = list(sub_st_meta_df.index)


    # mapping cells to spots, celltype by celltype
    mapped_cells, mapped_spots = [], []

    for tmp_ct in simulate_target_cts:
        ct_sc_cells = list(sub_sc_meta_df.loc[sub_sc_meta_df.celltype == tmp_ct,:].index)
        ct_st_spots = list(sub_st_meta_df.loc[sub_st_meta_df.celltype == tmp_ct,:].index)

        if len(ct_sc_cells) < len(ct_st_spots)*min_cell_per_spot:
            # in case use all cells of this celltype in sc data, but still not enough for the min cell per spot
            ct_sc_cells = ct_sc_cells * (len(ct_st_spots)*min_cell_per_spot//len(ct_sc_cells) + 1)

        # make sure the min cell number in a spot
        ct_mapped_cells = random.sample(ct_sc_cells, k=len(ct_st_spots)*min_cell_per_spot)
        ct_mapped_spots = random.sample(ct_st_spots*min_cell_per_spot, k=len(ct_st_spots)*min_cell_per_spot)

        # then randomly map the rest cells
        ct_rest_cells = [c for c in ct_sc_cells if c not in ct_mapped_cells]
        tmp_mapping_spots = ct_st_spots * (max_cell_per_spot - min_cell_per_spot)
        random.shuffle(tmp_mapping_spots)
        rest_mapped_spots = random.sample(tmp_mapping_spots, k=len(ct_rest_cells))

        mapped_cells += ct_mapped_cells + ct_rest_cells
        mapped_spots += ct_mapped_spots + rest_mapped_spots

    print('> {} cells are mapped'.format(len(mapped_cells)))


    simu_spot_cell_dic = {}
    for i, tmp_cell in enumerate(mapped_cells):
        if mapped_spots[i] not in simu_spot_cell_dic:
            simu_spot_cell_dic.update({
                mapped_spots[i]: [tmp_cell]
            })
        else:
            simu_spot_cell_dic[mapped_spots[i]].append(tmp_cell)

    # fill the simulated st count data
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


    
    ##########################
    # output simulation data #
    ##########################




    simulate_output_dir = simulate_output_dir
    
    # sc
    sc_count_df_enriched.to_csv('{}/sc_count_round1.tsv'.format(simulate_output_dir), sep='\t')
    sub_sc_meta_df.to_csv('{}/sc_meta.tsv'.format(simulate_output_dir), sep ='\t', header = None)

    # st
    simulated_st_count_df.to_csv('{}/st_count_round1.tsv'.format(simulate_output_dir), sep ='\t')
    sub_st_pos_df.to_csv('{}/st_coord.tsv'.format(simulate_output_dir), sep ='\t')
    sub_st_meta_df.to_csv('{}/st_meta.tsv'.format(simulate_output_dir), sep ='\t', header = None)

    # target ips 
    with open('{}/target_ct_pair_ip_dic.pkl'.format(simulate_output_dir),'wb') as f:
        pkl.dump(target_ct_pair_ip_dic, f)
        
        
    # ct_pair-ip-TF-target_gene dic
    with open('{}/target_ct_pair_ip_tf_tg_dic.pkl'.format(simulate_output_dir),'wb') as f:
        pkl.dump(target_ct_pair_ip_tf_tg_dic, f)

    # target celltype pair types
    tmp_output_df = pd.DataFrame(['near']*len(near_target_ct_pairs)+['far']*len(far_target_ct_pairs),
                                index = all_target_ct_pairs,
                                columns = ['type'])
    tmp_output_df.to_csv('{}/target_ct_pair_type.tsv'.format(simulate_output_dir), sep='\t')

    # simulate parameters
    tmp_output_df = pd.DataFrame([add_sigma, add_fc, nip_per_ctp, target_ip_cell_per, min_cells, min_cell_per_spot, max_cell_per_spot],
                                index = ['add_sigma', 'add_fc', 'nip_per_ctp', 'target_ip_cell_per', 'min_cells', 'min_cell_per_spot', 'max_cell_per_spot'],
                                columns = ['value'])
    tmp_output_df.to_csv('{}/simulate_parameters.tsv'.format(simulate_output_dir), sep='\t')

    # prepare celltype_frac.txt for stlearn input
    # we should regenerate celltype_frac.txt based on simulation mapping result
    target_ct_frac = pd.DataFrame(
        np.zeros((simulated_st_count_df.shape[1], len(simulate_target_cts))),
        index = list(simulated_st_count_df.columns),
        columns = simulate_target_cts
    )

    for tmp_spot in simu_spot_cell_dic:
        tmp_ct_sr = sub_sc_meta_df.loc[simu_spot_cell_dic[tmp_spot],:]
        tmp_frac_sr = tmp_ct_sr.value_counts()/tmp_ct_sr.shape[0]
        for ct in list(tmp_frac_sr.index):
            target_ct_frac.loc[tmp_spot,ct] = tmp_frac_sr.loc[ct,]
    target_ct_frac.to_csv('{}/st_spot_celltype_frac.txt'.format(simulate_output_dir), sep='\t')


    # output spot-cell mapping relation
    with open('{}/simu_spot_cell_dic.pkl'.format(simulate_output_dir),'wb') as f:
        pkl.dump(simu_spot_cell_dic, f)


    # output only target celltype pair distype
    select_index = [ctp in all_target_ct_pairs for ctp in list(ct_distype_sr.index)]
    target_ct_distype_sr = ct_distype_sr.loc[select_index]
    with open('{}/target_ct_distype_sr.pkl'.format(simulate_output_dir), 'wb') as f:
        pkl.dump(target_ct_distype_sr, f)

    # output all target ips
    def ip_name_trans(ip):
        ip = ip.split('_')
        if len(ip) == 2:
            return ' - '.join(ip)
        elif len(ip) > 2:
            return '{} - ({})'.format(ip[0], '+'.join(ip[1:]))
        else:
            print('ip: {} is not a standard ip'.format(ip))
            return None

    all_target_ips_name2 = []
    for ip_i, tmp_ip in enumerate(all_target_ips):
        tmp_ip_name2 = ip_name_trans(tmp_ip)
        if tmp_ip_name2 != None:
            all_target_ips_name2.append(tmp_ip_name2)
        else:
            all_target_ips.remove(tmp_ip)

    all_target_ips_df = pd.DataFrame({
        'v1': all_target_ips,
        'v2': all_target_ips_name2
    })
    all_target_ips_df.to_csv('{}/all_target_ips.tsv'.format(simulate_output_dir), 
                             sep='\t', index=None, header=None)

