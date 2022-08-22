import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import pickle as pkl
import os
import sys
from sklearn.metrics import roc_curve, auc, precision_recall_curve, f1_score
from sklearn.cluster import KMeans
import ot


def cal_ct_avg_dis(meta_df, pos_df, ct_name_correct_dic=None, nspot_cutoff=5):
    
    avg_dis_sr, ct_name_list = [], []
        
    if ct_name_correct_dic != None:
        meta_df['celltype'] = meta_df['celltype'].apply(lambda x: ct_name_correct_dic[x] if x in ct_name_correct_dic else x)
        
    ct_list = list(set(list(meta_df['celltype'])))  
    
    for n, ct_x in enumerate(ct_list):
        ct_x_barcode = list(meta_df.loc[meta_df.celltype == ct_x,:].index)
        
        if len(ct_x_barcode) < nspot_cutoff:
            continue
            
        select_barcode = [bar in ct_x_barcode for bar in list(pos_df.index)]
        sub_pos_df_x = pos_df.loc[select_barcode,:]

        for m in range(n+1, len(ct_list)):
            ct_y = ct_list[m]
            ct_y_barcode = list(meta_df.loc[meta_df.celltype == ct_y,:].index)
            
            if len(ct_y_barcode) < nspot_cutoff:
                continue
            
            select_barcode = [bar in ct_y_barcode for bar in list(pos_df.index)]
            sub_pos_df_y = pos_df.loc[select_barcode,:]
            
            ct_name_list.append('{}|{}'.format(ct_x, ct_y))
            
            tmp_dis_matrix = ot.dist(np.array(sub_pos_df_x.loc[:,['row','col']]), 
                                     np.array(sub_pos_df_y.loc[:,['row','col']]),
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



def plot_st_celltype(tmp_pos_df, tmp_meta_df,
                    plot_title = 'ST celltype',
                    meta_col = 'celltype',
                    pos_x_col = 'row',
                    pos_y_col = 'col',
                     s=20,
                    target_ct = None):

    tmp_ct_list = set(list(tmp_meta_df[meta_col]))
    color_list = sns.hls_palette(len(tmp_ct_list), l=0.45, s=0.85)
    
    if target_ct == None:
        target_ct = tmp_ct_list
    elif type(target_ct) == str:
        target_ct = [target_ct]
    elif type(target_ct) == list:
        target_ct = target_ct
    
    plt.figure(figsize=[5,5])

    for i, celltype in enumerate(sorted(tmp_ct_list)):
        select_index = [tmp_ct == celltype for tmp_ct in list(tmp_meta_df.loc[:,meta_col])]
        select_barcode = list(tmp_meta_df.loc[select_index,:].index)
        x_list = [int(tmp_pos_df.loc[barcode,pos_x_col]) for barcode in select_barcode]
        y_list = [int(tmp_pos_df.loc[barcode,pos_y_col]) for barcode in select_barcode]
        if celltype in target_ct:  
            plt.scatter(x_list,y_list,color=color_list[i], s=s, alpha=0.8)
        else:
            plt.scatter(x_list,y_list,color='gray', s=s, alpha=0.15)

    plt.xticks([])
    plt.yticks([])
    lg = plt.legend(sorted(tmp_ct_list),bbox_to_anchor=(1.0,0.75), frameon=False, fontsize=12)
    for n in range(len(tmp_ct_list)):
        lg.legendHandles[n]._sizes = [20]

    plt.title(plot_title,fontsize=16, pad=15, fontweight='bold')
    plt.show()


def extract_ip_cc(pval_path, ct_list, name_trans_path, ct_name_correct=None, cc_ip_path = None, ntop=None):
    
    if not os.path.exists(pval_path):
        print('cellchat result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic

    
    res_dic = {}
    ## key: ct;
    ## value: [ip_pair(str)]
    
    ip_pval_df = pd.read_csv(pval_path, sep = '\t')
    
    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])
    
    if ct_name_correct != None:        
        tmp_col = list(ip_pval_df.columns)
        for raw_name in ct_name_correct:
            tmp_col = [ct.replace(raw_name, ct_name_correct[raw_name]) for ct in tmp_col]
        ip_pval_df.columns = tmp_col
        
    ip_name_trans_dic = {}
    with open(name_trans_path,'r') as f:
        for line in f.readlines()[1:]:
            ip_name_trans_dic.update({line.strip().split('\t')[0]:line.strip().split('\t')[1]})


        
    for ct_type in ct_list:
        if ntop != None:
            tmp_pval_df = pd.DataFrame({'ip':[],'pval':[]})
        tmp_ip_list = []
        for ct in ['{}|{}'.format(ct_type.split('|')[1],ct_type.split('|')[0]), ct_type]:
            if ct not in list(ip_pval_df.columns):
                continue
            ct_pval_df = ip_pval_df.loc[:,ct]
            ct_pval_df = ct_pval_df.loc[ct_pval_df < 0.05]
            ### ct_pval_df to store the significant ips.
            
            if ntop != None:
                tmp_pval_df = tmp_pval_df.append(pd.DataFrame({'ip':list(ct_pval_df.index),'pval':list(ct_pval_df)}))
                
                
            sig_ips = list(ct_pval_df.index)
                
            sig_ips = [ip_name_trans_dic[ip] for ip in sig_ips]            
            tmp_ip_list += sig_ips
            
        if ntop != None:
            tmp_pval_df = tmp_pval_df.drop_duplicates()
            tmp_pval_df = tmp_pval_df.sort_values(by='pval',ascending=True)
            tmp_ip_list = [ip_name_trans_dic[ip] for ip in list(tmp_pval_df['ip'])]
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            tmp_ip_list = tmp_ip_list[:ntop]
        else:
            tmp_ip_list = list(set(tmp_ip_list))       
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            
        res_dic.update({ct_type:tmp_ip_list})
                
    return res_dic




def extract_ip_cpdb(pval_path, ct_list, ip_trans_path, dec_path, ct_name_correct = None, cc_ip_path = None, ntop=None):
    
    if not os.path.exists(pval_path):
        print('cellphonedb result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic

    
    res_dic = {}
    ## key: ct;
    ## value: [ip_pair(str)]

    ip_pval_df = pd.read_csv(pval_path,sep= '\t',index_col = 1)
    
    ip_name_trans_dic = {}
    with open(ip_trans_path,'r') as f:
        for line in f.readlines()[1:]:
            ip_name_trans_dic.update({line.strip().split('\t')[0]:line.strip().split('\t')[1]})

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])

    if ct_name_correct != None:
        tmp_col = list(ip_pval_df.columns)
        for raw_name in ct_name_correct:
            tmp_col = [ct.replace(raw_name, ct_name_correct[raw_name]) for ct in tmp_col]
        ip_pval_df.columns = tmp_col
        
    for ct_type in ct_list:
        tmp_ip_list = []
        if ntop != None:
            tmp_pval_df = pd.DataFrame({'ip':[],'pval':[]})

        for ct in ['{}|{}'.format(ct_type.split('|')[1],ct_type.split('|')[0]), ct_type]:
            if ct not in list(ip_pval_df.columns):
                continue
        
            ct_pval_df = ip_pval_df.loc[:,ct]
            ct_pval_df = ct_pval_df.loc[ct_pval_df < 0.05]                

            cpdb_name_trans_dic = {}
            with open(dec_path,'r') as f:
                for line in f.readlines()[1:]:
                    cpdb_name_trans_dic.update({line.strip().split('\t')[1]:line.strip().split('\t')[0]})
            new_cpdb_index = []
            for ip in list(ct_pval_df.index):
                tmp_list = []
                for gene in ip.split('_'):
                    if gene in cpdb_name_trans_dic:
                        gene = cpdb_name_trans_dic[gene]
                    tmp_list.append(gene)
                new_cpdb_index.append('_'.join(tmp_list))
            ct_pval_df.index = new_cpdb_index
            
            if ntop != None:
                tmp_pval_df = tmp_pval_df.append(pd.DataFrame({'ip':list(ct_pval_df.index),'pval':list(ct_pval_df)}))
            
            tmp_ip_list += [ip_name_trans_dic[ip] for ip in list(ct_pval_df.index) if ip in ip_name_trans_dic]
#         tmp_ip_list = list(set(tmp_ip_list))
        
        if ntop != None:
            tmp_pval_df = tmp_pval_df.drop_duplicates()
            tmp_pval_df = tmp_pval_df.sort_values(by='pval',ascending=True)
            tmp_ip_list = [ip_name_trans_dic[ip] for ip in list(tmp_pval_df['ip']) if ip in ip_name_trans_dic]
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            tmp_ip_list = tmp_ip_list[:ntop]
        else:
            tmp_ip_list = list(set(tmp_ip_list))       
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

#         if cc_ip_path != None:
#             tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
#         if ntop != None:
#             tmp_ip_list = tmp_ip_list[:ntop]

        res_dic.update({ct_type:tmp_ip_list})
    
    return res_dic


def extract_ip_italk(res_path,ct_list,ct_name_correct=None, cc_ip_path = None, ntop=None):

    if not os.path.exists(res_path):
        print('italk result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic

    res_dic = {}
    
    res_df = pd.read_csv(res_path,sep='\t')
    ## filter out 0 exp
    res_df = res_df.loc[res_df.cell_from_mean_exprs > 0,:]
    res_df = res_df.loc[res_df.cell_to_mean_exprs > 0, :]

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])

    if ct_name_correct != None:
        res_df['cell_from'] = res_df['cell_from'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)
        res_df['cell_to'] = res_df['cell_to'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)
    
    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        if ct_a in set(res_df['cell_from']) and ct_b in set(res_df['cell_to']):
            sub_df = res_df.loc[res_df.cell_from == ct_a,:].loc[res_df.cell_to == ct_b,:]
            sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis = 1)
            tmp_ip_list += list(sub_df['ip'])
        if ct_b in set(res_df['cell_from']) and ct_a in set(res_df['cell_to']):
            sub_df = res_df.loc[res_df.cell_from == ct_b,:].loc[res_df.cell_to == ct_a,:]
            sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis = 1)
            tmp_ip_list += list(sub_df['ip'])
            
        if ntop != None:
            if ct_a in set(res_df['cell_from']) and ct_b in set(res_df['cell_to']) and ct_b in set(res_df['cell_from']) and ct_a in set(res_df['cell_to']):
                tmp_score_df = res_df.loc[res_df.cell_from == ct_b,:].loc[res_df.cell_to == ct_a,:]
                tmp_score_df = tmp_score_df.append(res_df.loc[res_df.cell_from == ct_a,:].loc[res_df.cell_to == ct_b,:])
                tmp_score_df['ip'] = tmp_score_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis = 1)
                tmp_score_df['score'] = tmp_score_df['cell_from_mean_exprs'] * tmp_score_df['cell_to_mean_exprs']
                tmp_score_df = tmp_score_df.sort_values(by='score',ascending=False)
                tmp_ip_list = list(tmp_score_df['ip'])
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                tmp_ip_list = tmp_ip_list[:ntop]
        else:         
            tmp_ip_list = list(set(tmp_ip_list))   
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})
                
    return res_dic


def extract_ip_scr(res_dir,ct_list, lr_score=0.4, cc_ip_path = None, ntop=None):

    if not os.path.exists(res_dir):
        print('scr result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic

    res_dic = {}
    res_ct_list = [f.split('.')[0] for f in os.listdir(res_dir)]
    ## the output result file name list
    ## result file name format : Ductal_Endothelial.tsv

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    
    
    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        # ct_a for sender celltype
        # ct_b for recevier celltype
        # ct_a -> ct_b
        
        # read scr result file named by celltype pair            
        if '{}_{}'.format(ct_a, ct_b) in res_ct_list and '{}_{}'.format(ct_b, ct_a) in res_ct_list:
            ### confirm the result matrix has the specific celltype interaction type
            for pair_index, ct_pair in enumerate(['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]):

                        ####### extract predicting score. from tools result #######'
                sub_df = pd.read_csv('{}/{}_{}.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]),sep='\t')
                sub_df = sub_df.loc[sub_df.LRscore > lr_score,:]

                ## the scores of italk result are always larger than 0
                sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x[0],x[1]), axis = 1)
                
                if ntop != None:
                    if pair_index == 0:
                        tmp_score_df = sub_df
                    else:
                        tmp_score_df = tmp_score_df.append(sub_df)

                tmp_ip_list += list(sub_df['ip'])
                
            if ntop != None:
                tmp_score_df = tmp_score_df.drop_duplicates()
                tmp_score_df = tmp_score_df.sort_values(by='LRscore',ascending=False)
                tmp_ip_list = list(tmp_score_df['ip'])
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                tmp_ip_list = tmp_ip_list[:ntop]
            else:
                tmp_ip_list = list(set(tmp_ip_list))        
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})
                
    return res_dic


def extract_ip_cytotalk(res_dir,ct_list, cc_ip_path = None, ntop=None):
    
    if not os.path.exists(res_dir):
        print('cytotalk result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic

    
    res_dic = {}
    res_ct_list = os.listdir(res_dir)

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    
    
    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        
        if ntop != None:
            tmp_score_df = pd.DataFrame({})

        if '{}_{}'.format(ct_a, ct_b) in res_ct_list or '{}_{}'.format(ct_b, ct_a) in res_ct_list:
            ### confirm the result matrix has the specific celltype interaction type
            for ct_index, ct_pair in enumerate(['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]):
                ct_pair2 = ct_pair.replace('|','_')
                if ct_pair2 in res_ct_list:
                    file_ct = ct_pair2
                    try: 
                        sub_df = pd.read_csv('{}/{}/FinalNetwork.txt'.format(res_dir, file_ct),sep='\t')
                        sub_df = sub_df.loc[sub_df.is_ct_edge == True,:]
                        sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x['node1'].upper(), x['node2'].upper()), axis=1)
                        sub_df['ct_pair'] = sub_df.apply(lambda x: '{}|{}'.format(x['node1_type'], x['node2_type']), axis=1)
                        
                        #########################################
                        # codes for old version CytoTalk output #
                        #########################################
                        # sub_df = pd.read_csv('{}/{}/IllustratePCSF/PCSF_CrosstalkEdge.txt'.format(res_dir, file_ct),sep='\t')
                        # sub_df = sub_df.loc[sub_df.CrosstalkScore > 0,:]
                        #sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x[0].split('(')[0].strip(), x[1].split('(')[0].strip()), axis = 1)

                        if ntop != None:
                            tmp_score_df = tmp_score_df.append(sub_df)
                        tmp_ip_list += list(sub_df['ip'])
                    except:
                        continue
        
        if ntop != None and len(tmp_ip_list) > 0:
            tmp_score_df = tmp_score_df.drop_duplicates()
            tmp_score_df = tmp_score_df.sort_values(by='cost',ascending=False)
            tmp_ip_list = list(tmp_score_df['ip'])
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            tmp_ip_list = tmp_ip_list[:ntop]
        else:
            tmp_ip_list = list(set(tmp_ip_list))        
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})
                    
    return res_dic


def extract_ip_natmi(res_path,ct_list, ct_name_correct = None, cc_ip_path = None, ntop=None, weight_threshold = 0.1):
    
    if not os.path.exists(res_path):
        print('natmi result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic

    
    res_dic = {}
    res_df = pd.read_csv(res_path)
    res_df = res_df.loc[:,['Sending cluster','Ligand symbol','Receptor symbol','Target cluster','Edge average expression derived specificity']]
    res_df.columns = ['ct_from','ligand','receptor','ct_to','weight']
    res_df = res_df.loc[res_df.weight >= weight_threshold,:]

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    if ct_name_correct != None:
        res_df['ct_from'] = res_df['ct_from'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)
        res_df['ct_to'] = res_df['ct_to'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)

    
    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        # ct_a for sender celltype
        # ct_b for recevier celltype
        # ct_a -> ct_b

        # read scr result file named by celltype pair            
        if ct_a in set(res_df.ct_from) and ct_b in set(res_df.ct_to):
            ### confirm the result matrix has the specific celltype interaction type
            for ct_index, ct_pair in enumerate(['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]):
                sub_df = res_df.loc[res_df.ct_from == ct_pair.split('|')[0],:].loc[res_df.ct_to == ct_pair.split('|')[1],:]
                if sub_df.shape[0] == 0:
                    continue

                sub_df['ip'] = sub_df.apply(lambda x: '{} - {}'.format(x[1], x[2]), axis = 1)
                
                if ntop != None:
                    if ct_index == 0:
                        tmp_score_df = sub_df
                    else:
                        tmp_score_df = tmp_score_df.append(sub_df)
                
                tmp_ip_list += list(sub_df['ip'])
        
        if ntop != None and ct_a in set(res_df.ct_from) and ct_b in set(res_df.ct_to):
            tmp_score_df = tmp_score_df.drop_duplicates()
            tmp_score_df = tmp_score_df.sort_values(by='weight',ascending=False)
            tmp_ip_list = list(tmp_score_df['ip'])
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            tmp_ip_list = tmp_ip_list[:ntop]
        else:            
            tmp_ip_list = list(set(tmp_ip_list))        
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
        
        res_dic.update({ct_type:tmp_ip_list})

                
    return res_dic


def trans_symbol(ip):
    ## be used in extract_ip_icellnet
    if ip.count('+') > 0:
        tmp_ip = ip.split(' - ')
        ip = ' - '.join([tmp_ip[0], '({})'.format('+'.join(tmp_ip[-1].split(' + ')))])
    return ip

def extract_ip_icellnet(res_dir,ct_list, lr_score=20, cc_ip_path = None, ntop=None):
    ### ct_distype & ip_dis_df are calculated in PAAD data prepare part
    res_dic = {} 
    res_ct_list = [f.split('.')[0][6:] for f in os.listdir(res_dir)]

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        
        if ntop != None:
            tmp_score_df = pd.DataFrame({})
            # split df by celltype
        if '{}_{}'.format(ct_a,ct_b) in res_ct_list or '{}_{}'.format(ct_b,ct_a) in res_ct_list:
            ### confirm the result matrix has the specific celltype interaction type
            for ct_index, ct_pair in enumerate(['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]):
                if ct_pair in ct_list:
                    sub_df = pd.read_csv('{}/score_{}_{}.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]), sep='\t', header=None)
                    sub_df.columns = ['ip','score']
                    sub_df.dropna(how='any', inplace=True)
                    sub_df = sub_df.loc[sub_df.score > lr_score,:]

                    sub_df['ip'] = sub_df['ip'].apply(lambda x: x.replace(' / ',' - '))
                    sub_df['ip'] = sub_df['ip'].apply(trans_symbol)
                    
                    if ntop != None:
                        tmp_score_df = tmp_score_df.append(sub_df)

                    tmp_ip_list += list(sub_df['ip'])
        
        if ntop != None and ('{}_{}'.format(ct_a,ct_b) in res_ct_list or '{}_{}'.format(ct_b,ct_a) in res_ct_list):
            tmp_score_df = tmp_score_df.drop_duplicates()
            if tmp_score_df.shape[0] == 0:
                tmp_ip_list = []
            else:
                tmp_score_df = tmp_score_df.sort_values(by='score',ascending=False)
                tmp_ip_list = list(tmp_score_df['ip'])
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                tmp_ip_list = tmp_ip_list[:ntop]
        else:            
            tmp_ip_list = list(set(tmp_ip_list))        
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})

    return res_dic
  
def extract_ip_icellnet_refine(res_dir,ct_list, cc_ip_path = None, lr_score_cutoff = 20,ntop=None,pair_top=100):
    
    if not os.path.exists(res_dir):
        print('icellnet result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic

    
    res_dic = {} 
    res_ct_file_list = [f for f in os.listdir(res_dir) if f.startswith('score') and f.count('pdf')==0]
    res_ct_list = [f.replace('score_','').replace('_out_all.tsv','') for f in res_ct_file_list]
    

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    for fn, ct_type in enumerate(res_ct_list):
        tmp_ip_list = []
        ct_a = ct_type
        
        pair_df = pd.read_csv('{}/score_{}_out_all.tsv'.format(res_dir,ct_a),sep='\t')
        pair_df = pair_df.iloc[:pair_top,:]
        
        ct_b_list = list(pair_df.columns)
        
        for ct_b in ct_b_list:
                
            if ntop != None:
                tmp_score_df = pd.DataFrame({})
                # split df by celltype
            if '{}|{}'.format(ct_a,ct_b) in ct_list:
                ### confirm the result matrix has the specific celltype interaction type
                ct_pair = '{}|{}'.format(ct_a,ct_b)
                if ct_pair in ct_list:
                    sub_df = pair_df.loc[:,ct_b]
                    sub_df.fillna(0,inplace=True)
                    sub_df = sub_df.loc[sub_df>lr_score_cutoff]
                    sub_ip_list = list(sub_df.index)
                    sub_ip_list = [x.replace(' / ',' - ') for x in sub_ip_list]
                    sub_ip_list = [trans_symbol(x) for x in sub_ip_list]
#                         sub_df = pd.read_csv('{}/score_{}_{}.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]), sep='\t', header=None)
#                     sub_df.columns = ['ip','score']
#                     sub_df.dropna(how='any', inplace=True)
#                     sub_df = sub_df.loc[sub_df.score > lr_score,:]

#                     sub_df['ip'] = sub_df['ip'].apply(lambda x: x.replace(' / ',' - '))
#                     sub_df['ip'] = sub_df['ip'].apply(trans_symbol)

                    if ntop != None:
                        tmp_score_df = tmp_score_df.append(pd.DataFrame({'ip':sub_ip_list,'score':list(sub_df)}))

                    tmp_ip_list += sub_ip_list
        
            if ntop != None and '{}|{}'.format(ct_a,ct_b) in ct_list:
                tmp_score_df = tmp_score_df.drop_duplicates()
                if tmp_score_df.shape[0] == 0:
                    tmp_ip_list = []
                else:
                    tmp_score_df = tmp_score_df.sort_values(by='score',ascending=False)
                    tmp_ip_list = list(tmp_score_df['ip'])
                    if cc_ip_path != None:
                        tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                    tmp_ip_list = tmp_ip_list[:ntop]
            else:            
                tmp_ip_list = list(set(tmp_ip_list))        
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                    
            for ct_dic_type in ['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]:
                if ct_dic_type in ct_list:
                    if ct_dic_type not in res_dic:
                        res_dic.update({ct_dic_type:tmp_ip_list})
                    else:
                        res_dic[ct_dic_type] += tmp_ip_list

#             res_dic.update({ct_type:tmp_ip_list})
    for ct in ct_list:
        if ct not in res_dic:
            res_dic.update({ct:[]})

    return res_dic

def extract_ip_nichenet(res_dir,ct_list, cc_ip_path = None, ntop=None):
    
    if not os.path.exists(res_dir):
        print('nichenet result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic

    res_dic = {}
    res_ct_list = os.listdir(res_dir)

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        
        if ntop != None:
            tmp_score_list = []   
            
        if 'output{}_{}'.format(ct_a, ct_b) in res_ct_list:
            ### confirm the result matrix has the specific celltype interaction type
            for ct_pair in ['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]:
                if ct_pair in ct_list:
                    sub_df = pd.read_csv('{}/output{}_{}/LR_potential.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]),sep='\t')
                    for l_i in range(sub_df.shape[0]):
                        for c_j in range(sub_df.shape[1]):
                            if sub_df.iloc[l_i,c_j] != 0:
                                tmp_ip_list.append('{} - {}'.format(sub_df.columns[c_j], sub_df.index[l_i]))
                                if ntop != None:
                                    tmp_score_list.append(sub_df.iloc[l_i,c_j])
        
        if ntop != None:
            tmp_score_sr = pd.Series(tmp_score_list, index= tmp_ip_list)
            tmp_score_sr = tmp_score_sr.sort_values(ascending=False)
            tmp_ip_list = list(tmp_score_sr.index)
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
            tmp_ip_list = tmp_ip_list[:ntop]
        else:
            tmp_ip_list = list(set(tmp_ip_list))        
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})

                
    return res_dic
    

def extract_ip_scmlnet(res_dir,ct_list, cc_ip_path = None, ntop=None):
    
    if not os.path.exists(res_dir):
        print('scmlnet result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic
    
    res_dic = {}
    res_ct_list = os.listdir(res_dir)

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        
#         if ntop != None:
#             tmp_score_df = pd.DataFrame({})
        
            # split df by celltype
        for ct_pair in ['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]:
            if '{}_{}_RecClus.tsv'.format(ct_pair.split('|')[0], ct_pair.split('|')[1]) in res_ct_list and '{}_{}_LigClus.tsv'.format(ct_pair.split('|')[0], ct_pair.split('|')[1]) in res_ct_list and '{}_{}_LR.txt'.format(ct_pair.split('|')[0], ct_pair.split('|')[1]) in res_ct_list:
            ### confirm the result matrix has the specific celltype interaction type
                with open('{}/{}_{}_LR.txt'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]), 'r') as f:
                    lr_list = [line.strip() for line in f.readlines()]
                tmp_r_df = pd.read_csv('{}/{}_{}_RecClus.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]), sep='\t')
                tmp_l_df = pd.read_csv('{}/{}_{}_LigClus.tsv'.format(res_dir, ct_pair.split('|')[0], ct_pair.split('|')[1]), sep='\t')
                tmp_ip_list += [lr.replace('_',' - ') for lr in lr_list]
        tmp_ip_list = list(set(tmp_ip_list))
        
        if cc_ip_path != None:
            tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
        if ntop != None:
            tmp_ip_list = tmp_ip_list[:ntop]

        res_dic.update({ct_type:tmp_ip_list})
                
    return res_dic


def extract_ip_connectome(res_path, ct_list, ct_name_correct = None, cc_ip_path = None, ntop=None):
    
    if not os.path.exists(res_path):
        print('connectome result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic

    
    res_dic = {}
    res_df = pd.read_csv(res_path,sep='\t',header=None)
    res_df.columns = ['edge','score']
    res_df['ct_from'] = res_df['edge'].apply(lambda x: x.split(' - ')[0])
    res_df['ct_to'] = res_df['edge'].apply(lambda x: x.split(' - ')[-1])
    res_df['ip'] = res_df['edge'].apply(lambda x: '{} - {}'.format(x.split(' - ')[1], x.split(' - ')[2]))
    ## filter out 0 exp
    res_df = res_df.loc[res_df.score > 0,:]

    if cc_ip_path != None:
        cc_ip_df = pd.read_csv(cc_ip_path, sep='\t')
        ip_list_all = list(cc_ip_df['interaction_name_2'])    

    if ct_name_correct != None:
        res_df['ct_from'] = res_df['ct_from'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)
        res_df['ct_to'] = res_df['ct_to'].apply(lambda x: ct_name_correct[x] if x in ct_name_correct else x)
    
    for ct_type in ct_list:
        tmp_ip_list = []
        ct_a, ct_b = ct_type.split('|')
        
        if ntop != None:
            tmp_score_df = pd.DataFrame({})
        
        if ct_a in set(res_df['ct_from']) and ct_b in set(res_df['ct_to']):
            ### confirm the result matrix has the specific celltype interaction type
            for ct_pair in ['{}|{}'.format(ct_a,ct_b), '{}|{}'.format(ct_b,ct_a)]:
                if ct_pair in ct_list:

                    sub_df = res_df.loc[res_df.ct_from == ct_pair.split('|')[0],:].loc[res_df.ct_to == ct_pair.split('|')[1],:]
                    if sub_df.shape[0] == 0:
                        continue
                    
                    if ntop != None:
                        tmp_score_df = tmp_score_df.append(sub_df)
                        
                    tmp_ip_list += list(sub_df['ip'])
                    
        
        if ntop != None:
            if tmp_score_df.shape[0] == 0:
                tmp_ip_list = []
            else:
                tmp_score_df = tmp_score_df.drop_duplicates()
                tmp_score_df = tmp_score_df.sort_values(by='score',ascending=False)
                tmp_ip_list = list(tmp_score_df['ip'])
                if cc_ip_path != None:
                    tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]
                tmp_ip_list = tmp_ip_list[:ntop]
        else:            
            tmp_ip_list = list(set(tmp_ip_list))        
            if cc_ip_path != None:
                tmp_ip_list = [ip for ip in tmp_ip_list if ip in ip_list_all]

        res_dic.update({ct_type:tmp_ip_list})
                
    return res_dic


def extract_ip_cellcall(res_path_cellcall, tmp_ct_list, score_threshold=0):
    
    if not os.path.exists(res_path_cellcall):
        print('cellcall result not available')
        tmp_res_dic = {}
        for tmp_ctp in tmp_ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic
    
    tmp_res_df = pd.read_csv(res_path_cellcall, sep='\t', index_col = 0)
    tmp_target_ct_pairs = tmp_ct_list + ['|'.join(list(reversed(ctp.split('|')))) for ctp in tmp_ct_list]
    tmp_res_df = tmp_res_df.loc[:,tmp_target_ct_pairs]
    
    tmp_res_dic = {}
    for ctp_i, tmp_ctp in enumerate(tmp_ct_list):
        
        # subset result df by target celltypes (two directions)
        tmp_ctp_rev = tmp_target_ct_pairs[ctp_i+len(tmp_ct_list)]
        tmp_sub_res_df = tmp_res_df.loc[:,[tmp_ctp, tmp_ctp_rev]]
        
        # select interactions with score > score_threshold 
        tmp_sub_res_df = tmp_sub_res_df.loc[(tmp_sub_res_df > score_threshold).sum(axis=1) > 0,:]
        
        # check ip names
        tmp_ips = []
        for tmp_ip in list(tmp_sub_res_df.index):
            if tmp_ip.count(',') == 0: # single subunit
                tmp_ips.append(tmp_ip)
            else: # multi-subunits
                tmp_ip = tmp_ip.replace(',','+')
                tmp_ip = '{} - ({})'.format(tmp_ip.split(' - ')[0], tmp_ip.split(' - ')[1])
                tmp_ips.append(tmp_ip)
        
        # save results
        tmp_res_dic.update({
            tmp_ctp: tmp_ips
        })
    
    return tmp_res_dic


def extract_ip_domino(res_dir_domino, tmp_ct_list, avail_ip_list, score_threshold=0):
    
    if not os.path.exists(res_dir_domino):
        print('domino result not available')
        tmp_res_dic = {}
        for tmp_ctp in tmp_ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic
    
    tmp_raw_ip_dic = {}
    for sub_file in os.listdir(res_dir_domino):
        if sub_file.endswith('tsv'):
            tmp_sub_df = pd.read_csv('{}/{}'.format(res_dir_domino, sub_file), sep='\t', index_col = 0)
            for ctp_i, tmp_ctp in enumerate(list(tmp_sub_df.columns)):
                tmp_ips = list(tmp_sub_df.loc[tmp_sub_df.iloc[:,ctp_i] > score_threshold,:].index)
                tmp_raw_ip_dic.update({
                    tmp_ctp: tmp_ips
                })
    
    tmp_res_dic = {}
    for tmp_ctp in tmp_ct_list:
        tmp_ctp_rev = '|'.join(list(reversed(tmp_ctp.split('|'))))
        tmp_ips = []
        for ctp in [tmp_ctp, tmp_ctp_rev]:
            if ctp in tmp_raw_ip_dic:
                tmp_ips += tmp_raw_ip_dic[ctp]
        tmp_ips = list(set(tmp_ips).intersection(set(avail_ip_list)))
        tmp_res_dic.update({
            tmp_ctp: tmp_ips
        })
    
    return tmp_res_dic


def extract_ip_stlearn(res_path_stlearn, ct_list, score_threshold=0):
    
    if not os.path.exists(res_path_stlearn):
        print('stlearn result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic
    
    tmp_res_df = pd.read_csv(res_path_stlearn, sep='\t', index_col = 0)
    tmp_res_df.index = [tmp_ip.replace('_',' - ') for tmp_ip in list(tmp_res_df.index)]
    
    tmp_res_dic = {}
    for tmp_ctp in ct_list:
        tmp_ips = []
        for ctp in [tmp_ctp, '|'.join(list(reversed(tmp_ctp.split('|'))))]:
            if ctp in list(tmp_res_df.columns):
                tmp_ips += list(tmp_res_df.loc[tmp_res_df.loc[:,ctp] > score_threshold,:].index)
        tmp_ips = list(set(tmp_ips))
        tmp_res_dic.update({
            tmp_ctp: tmp_ips
        })
    
    return tmp_res_dic


def extract_ip_giotto(res_path_giotto, ct_list):
    
    if not os.path.exists(res_path_giotto):
        print('giotto result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic

    res_df = pd.read_csv(res_path_giotto, sep='\t')
    
    tmp_res_dic = {}
    for tmp_ctp in ct_list:
        tmp_ips = []
        for ctp in [tmp_ctp, '|'.join(list(reversed(tmp_ctp.split('|'))))]:
            ct_a, ct_b = ctp.split('|')
            sub_res_df = res_df.loc[res_df.lig_cell_type == ct_a,:].loc[res_df.rec_cell_type == ct_b,:]
            if sub_res_df.shape[0] != 0:                
                sub_res_df['ips'] = sub_res_df.apply(lambda x: '{} - {}'.format(x[3], x[6]), axis=1)
                tmp_ips += list(sub_res_df.loc[:,'ips'])
        tmp_ips = list(set(tmp_ips))
        tmp_res_dic.update({
            tmp_ctp: tmp_ips
        })
    
    return tmp_res_dic



def extract_ip_base_line(res_pkl_base_line, ct_list):
    
    if not os.path.exists(res_pkl_base_line):
        print('base line result not available')
        tmp_res_dic = {}
        for tmp_ctp in ct_list:
            tmp_res_dic.update({
                tmp_ctp: []
            })
        return tmp_res_dic
    
    tmp_res_dic = {}
    
    with open(res_pkl_base_line, 'rb') as f:
        tool_res_dic = pkl.load(f)
        
    for tmp_ctp in ct_list:
        tmp_ips = []
        for ctp in [tmp_ctp, '|'.join(list(reversed(tmp_ctp.split('|'))))]:
            if ctp in tool_res_dic:
                tmp_ips += tool_res_dic[ctp]
        tmp_ips = list(set(tmp_ips))
        tmp_res_dic.update({
            tmp_ctp: tmp_ips
        })
        
    return tmp_res_dic
    


def cal_es_score(rank_list,ip_set,rank_weight_list):
    
    ## generate the flag (whether the ips of ip_set in the rank_list)
    rank_weight_list = pd.Series(rank_weight_list)
    match_flag = pd.Series([1 if ip in ip_set else 0 for ip in rank_list])
    unmatch_flag = 1 - match_flag

    match_weight_sum = sum(rank_weight_list[match_flag == 1])
    if match_weight_sum != 0:
        match_weight_tag = 1/match_weight_sum
    else:
        match_weight_tag = 1

    unmatch_weight_sum = sum(rank_weight_list[unmatch_flag == 1])
    if unmatch_weight_sum != 0:
        unmatch_weight_tag = 1/unmatch_weight_sum
    else:
        unmatch_weight_tag = 1

    cumsum_score = np.cumsum(match_flag*rank_weight_list*match_weight_tag - unmatch_flag*rank_weight_list*unmatch_weight_tag)

    ES_score = max(cumsum_score) if max(cumsum_score) > -min(cumsum_score) else min(cumsum_score)
    
    return ES_score, cumsum_score


def cal_es_score_gsea(rank_list,ip_set,rank_weight_list):
    ## generate the flag (whether the ips of ip_set in the rank_list)
    rank_weight_list = pd.Series(rank_weight_list)
    match_flag = pd.Series([1 if ip in ip_set else 0 for ip in rank_list])
    unmatch_flag = 1 - match_flag

    match_weight_sum = sum(rank_weight_list[match_flag == 1])
    if match_weight_sum != 0:
        match_weight_tag = 1/match_weight_sum
    else:
        match_weight_tag = 1
        ### acutally it doesn't matter whether the match_weight_tag == 1, since the weight == 0

    unmatch_weight_sum = sum(rank_weight_list[unmatch_flag == 1])
    if sum(unmatch_flag) != 0:
        unmatch_weight_tag = 1/sum(unmatch_flag)
    else:
        unmatch_weight_tag = 1

    cumsum_score = np.cumsum(match_flag*rank_weight_list*match_weight_tag - unmatch_flag*unmatch_weight_tag)

    ES_score = max(cumsum_score) if max(cumsum_score) > -min(cumsum_score) else min(cumsum_score)
    
    return ES_score, cumsum_score


def generate_ES_dic_with_cumsum(d_rat_df, tool_res_dic, ct_distype_sr, 
                                origin=True, top_per=0.1, 
                                pval_cutoff_flag=False, drat_pval_cutoff=0.01):
    ES_dic, cumsum_dic = {},{}
    # {
    #     'tool':{
    #         'short':[ES_score],
    #         'long':[ES_score],
    #         'all':[ES_score]
    #     }
    # }

    # near/far celltype pairs
    near_ct_list = list(ct_distype_sr.loc[ct_distype_sr == 0].index)
    far_ct_list = list(ct_distype_sr.loc[ct_distype_sr == 1].index)
    
    # filtering short/long-range interactions
    max_index = int(d_rat_df.shape[0]*top_per)
    d_rat_df = d_rat_df.sort_values(by=['d_rat','p_val'])
    if pval_cutoff_flag:
        tmp_d_rat_df = d_rat_df.loc[d_rat_df.p_val <= drat_pval_cutoff,:]
        max_index = max(max_index, tmp_d_rat_df.shape[0])
    short_ip_df = d_rat_df.iloc[:max_index,:]

    # filtering long-range interactions
    d_rat_df = d_rat_df.sort_values(by=['d_rat','p_val'],ascending=[False,True])
    if pval_cutoff_flag:
        tmp_d_rat_df = d_rat_df.loc[d_rat_df.p_val <= drat_pval_cutoff,:]
        max_index = max(max_index, tmp_d_rat_df.shape[0])
    long_ip_df = d_rat_df.iloc[:max_index,:]


    tool_list = list(tool_res_dic.keys())
    
    
    for tool in tool_list:
        print(tool)
        tmp_short_es_list, tmp_long_es_list = [], []
        tmp_short_cs_list, tmp_long_cs_list = [], []
        tmp_short_ct_list, tmp_long_ct_list = [], []

        tmp_ct_list = list(tool_res_dic[tool].keys())

        for tmp_ctp in tool_res_dic[tool]:
            
            if tmp_ctp in near_ct_list:
                ## using short rank list
                if 'ip' in list(short_ip_df.columns):
                    s_rank_ip_list = list(short_ip_df['ip'])
                else:
                    s_rank_ip_list = list(short_ip_df.index)
                s_weight_list = list(1-short_ip_df['p_val'])
                
                if origin == True:
                    ES_score, cumsum_score = cal_es_score(
                        s_rank_ip_list,
                        tool_res_dic[tool][tmp_ctp],
                        s_weight_list)
                else:
                    ES_score, cumsum_score = cal_es_score_gsea(
                        s_rank_ip_list,
                        tool_res_dic[tool][tmp_ctp],
                        s_weight_list)
                tmp_short_es_list.append(ES_score)
                tmp_short_cs_list.append(cumsum_score)
                tmp_short_ct_list.append(tmp_ctp)

            elif tmp_ctp in far_ct_list:
                ## using long rank list
                if 'ip' in list(long_ip_df.columns):     
                    l_rank_ip_list = list(long_ip_df['ip'])
                else:
                    l_rank_ip_list = list(long_ip_df.index)
                l_weight_list = list(1-long_ip_df['p_val'])
                
                if origin == True:
                    ES_score, cumsum_score = cal_es_score(
                        l_rank_ip_list,
                        tool_res_dic[tool][tmp_ctp],
                        l_weight_list)
                else:
                    ES_score, cumsum_score = cal_es_score_gsea(
                        l_rank_ip_list, 
                        tool_res_dic[tool][tmp_ctp],
                        l_weight_list)
                tmp_long_es_list.append(ES_score)
                tmp_long_cs_list.append(cumsum_score)
                tmp_long_ct_list.append(tmp_ctp)

        ES_dic.update({tool:{
            'short':tmp_short_es_list,
            'short_ct':tmp_short_ct_list,
            'long':tmp_long_es_list,
            'long_ct':tmp_long_ct_list,
            'all':tmp_short_es_list+tmp_long_es_list,
            'all_ct':tmp_short_ct_list+tmp_long_ct_list
        }})

        cumsum_dic.update({tool:{
            'short':tmp_short_cs_list,
            'short_ct':tmp_short_ct_list,
            'long':tmp_long_cs_list,
            'long_ct':tmp_long_ct_list,
            'all':tmp_short_cs_list+tmp_long_cs_list,
            'all_ct':tmp_short_ct_list+tmp_long_ct_list
        }})

    return ES_dic, cumsum_dic



def generate_simulate_ES_dic_with_cumsum(d_rat_df, tool_res_dic, target_ct_pair_ip_dic,ct_distype_sr,
                                         top_per=0.1, origin_ls=False, origin_gsea=True,  
                                         pval_cutoff_flag=False, drat_pval_cutoff=0.01):
    ES_dic, cumsum_dic = {},{}
    
    # near/far celltype pairs
    near_ct_list = list(ct_distype_sr.loc[ct_distype_sr == 0].index)
    far_ct_list = list(ct_distype_sr.loc[ct_distype_sr == 1].index)

    # filtering short/long-range interactions if needed
    if origin_ls == True:
        # filtering short-range interactions
        max_index = int(d_rat_df.shape[0]*top_per)
        d_rat_df = d_rat_df.sort_values(by=['d_rat','p_val'])
        if pval_cutoff_flag:
            tmp_d_rat_df = d_rat_df.loc[d_rat_df.p_val <= drat_pval_cutoff,:]
            max_index = max(max_index, tmp_d_rat_df.shape[0])
        short_ip_df = d_rat_df.iloc[:max_index,:]

        # filtering long-range interactions
        d_rat_df = d_rat_df.sort_values(by=['d_rat','p_val'],ascending=[False,True])
        if pval_cutoff_flag:
            tmp_d_rat_df = d_rat_df.loc[d_rat_df.p_val <= drat_pval_cutoff,:]
            max_index = max(max_index, tmp_d_rat_df.shape[0])
        long_ip_df = d_rat_df.iloc[:max_index,:]
        
        
    # generate target ip rank for each celltype pair
    
    # if origin_ls (original long/short interactions) == True, 
    # we will add data-generated short/long interactions to simulated target interactions 
    target_ct_pair_ip_rank_dic = {}
    for tmp_ctp in target_ct_pair_ip_dic:
        tmp_target_ips = target_ct_pair_ip_dic[tmp_ctp]
        # change interaction names to interaction_name_2 type
        tmp_target_ips = []
        for tmp_ip in target_ct_pair_ip_dic[tmp_ctp]:
            tmp_ip = tmp_ip.split('_')
            if len(tmp_ip) == 2: # single subunit
                tmp_ip = ' - '.join(tmp_ip)
            else: # multi-subunits
                tmp_ip = '{} - ({})'.format(tmp_ip[0], '+'.join(tmp_ip[1:]))
            tmp_target_ips.append(tmp_ip)
        
        # rank target ips
        tmp_target_d_rat_df = d_rat_df.loc[tmp_target_ips,:]
        
        if origin_ls == True: # add short/long-range interactions
            if tmp_ctp in near_ct_list:
                tmp_target_d_rat_df = tmp_target_d_rat_df.append(short_ip_df)
                tmp_target_d_rat_df.drop_duplicates(inplace=True)
            elif tmp_ctp in far_ct_list:
                tmp_target_d_rat_df = tmp_target_d_rat_df.append(long_ip_df)
                tmp_target_d_rat_df.drop_duplicates(inplace=True)
            else:
                continue
        
        if tmp_ctp in near_ct_list:
            tmp_target_d_rat_df = tmp_target_d_rat_df.sort_values(by=['d_rat','p_val'])
        elif tmp_ctp in far_ct_list:
            tmp_target_d_rat_df = tmp_target_d_rat_df.sort_values(by=['d_rat','p_val'],ascending=[False,True])
        else:
            continue
        
        if 'ip' in list(tmp_target_d_rat_df.columns):
            tmp_target_ip_ranks = list(tmp_target_d_rat_df.loc[:,'ip'])
        else:
            tmp_target_ip_ranks = list(tmp_target_d_rat_df.index)
        
        target_ct_pair_ip_rank_dic.update({
            tmp_ctp:{
                'rank': tmp_target_ip_ranks,
                'weight': list(1-tmp_target_d_rat_df.loc[:,'p_val'])
            }
        })
    
    
    tool_list = list(tool_res_dic.keys())

    # calculate average ES for each tool
    for tool in tool_list:
        print(tool)
        tmp_short_es_list, tmp_long_es_list = [], [] # ES
        tmp_short_cs_list, tmp_long_cs_list = [], [] # cumsum ES
        tmp_short_ct_list, tmp_long_ct_list = [], [] # celltypes

        for tmp_ctp in tool_res_dic[tool]:
            if origin_gsea == True:
                ES_score, cumsum_score = cal_es_score(
                    target_ct_pair_ip_rank_dic[tmp_ctp]['rank'], 
                    tool_res_dic[tool][tmp_ctp] , 
                    target_ct_pair_ip_rank_dic[tmp_ctp]['weight'])
            else:
                ES_score, cumsum_score = cal_es_score_gsea(
                    target_ct_pair_ip_rank_dic[tmp_ctp]['rank'], 
                    tool_res_dic[tool][tmp_ctp] , 
                    target_ct_pair_ip_rank_dic[tmp_ctp]['weight'])
            
            if tmp_ctp in near_ct_list:
                tmp_short_es_list.append(ES_score)
                tmp_short_cs_list.append(cumsum_score)
                tmp_short_ct_list.append(tmp_ctp)
            elif tmp_ctp in far_ct_list:
                tmp_long_es_list.append(ES_score)
                tmp_long_cs_list.append(cumsum_score)
                tmp_long_ct_list.append(tmp_ctp)
                

        ES_dic.update({tool:{
            'short':tmp_short_es_list,
            'short_ct':tmp_short_ct_list,
            'long':tmp_long_es_list,
            'long_ct':tmp_long_ct_list,
            'all':tmp_short_es_list+tmp_long_es_list,
            'all_ct':tmp_short_ct_list+tmp_long_ct_list
        }})

        cumsum_dic.update({tool:{
            'short':tmp_short_cs_list,
            'short_ct':tmp_short_ct_list,
            'long':tmp_long_cs_list,
            'long_ct':tmp_long_ct_list,
            'all':tmp_short_cs_list+tmp_long_cs_list,
            'all_ct':tmp_short_ct_list+tmp_long_ct_list
        }})

    return ES_dic, cumsum_dic




def plot_ES_scatter(ES_dic, fig_save_dir, tool_name_trans_dic, 
                    h=6, w=2.5, s=30, marker='D', color='orange', 
                    fontsize=10, avg=False, figtail='', fontsize_x=12,
                    title_flag=True,
                    fig_save_flag=True, verbose=False, color_diff=False):
    
    tool_list = [tool_name_trans_dic[tool] for tool in list(ES_dic.keys())]
    
    short_avg_list, long_avg_list, all_avg_list = [],[],[]
    for tool in ES_dic:
        if len(ES_dic[tool]['short']) > 0 and len(ES_dic[tool]['long']) > 0:
            short_avg_list.append(np.average([ES_dic[tool]['short']]))
            long_avg_list.append(np.average([ES_dic[tool]['long']]))
            all_avg_list.append((np.average(ES_dic[tool]['short'])+np.average(ES_dic[tool]['long']))/2)
        elif len(ES_dic[tool]['short']) == 0:
            short_avg_list.append(-1)
            long_avg_list.append(np.average([ES_dic[tool]['long']]))
            all_avg_list.append(np.average(ES_dic[tool]['long']))
        elif len(ES_dic[tool]['long']) == 0:
            short_avg_list.append(np.average([ES_dic[tool]['short']]))
            long_avg_list.append(-1)
            all_avg_list.append(np.average(ES_dic[tool]['short']))
            
    short_avg_list = pd.Series(short_avg_list, index=tool_list)
    long_avg_list = pd.Series(long_avg_list, index=tool_list)
    all_avg_list = pd.Series(all_avg_list, index=tool_list)
        
    avg_rank = all_avg_list.rank(ascending=False)
    avg_rank = avg_rank.sort_values(ascending=False)
    rank_list = list(avg_rank.index)
    
    if fig_save_flag:
        if not os.path.exists(fig_save_dir):
            os.makedirs(fig_save_dir)

    plt.figure(figsize=(w,h))
    
    if avg == False:
        plt.scatter(short_avg_list.loc[rank_list], range(len(short_avg_list)), color='deepskyblue', marker='o', s=s, label='Short-range')
        plt.scatter(long_avg_list.loc[rank_list], range(len(long_avg_list)), color='orange', marker='D', s=s, label='Long-range')
        plt.yticks(range(len(long_avg_list)), rank_list, fontsize=fontsize)
        
        plt.legend(fontsize=8)
        plt.title('Distance Enrichment Score', fontsize=fontsize+2, pad=10, fontweight='bold')
        if fig_save_flag == True:
            plt.savefig('{}/tool_ES{}.pdf'.format(fig_save_dir,figtail), pad_inches = 0.1, bbox_inches = 'tight', dpi=500)
        plt.show()
    else:
        if color_diff:
            color_list = ['orange' if round(tmp_es,6) > -1 else 'gray' for tmp_es in list(all_avg_list.loc[rank_list])]
        else:
            color_list = ['orange']*len(all_avg_list)
        plt.scatter(all_avg_list.loc[rank_list], range(len(all_avg_list)), color=color_list, marker='D', s=s)
        plt.yticks(range(len(long_avg_list)), rank_list, fontsize=fontsize)
        if title_flag:
            plt.title('DES', fontsize=fontsize+2, pad=10, fontweight='bold')
        plt.xticks(fontsize=fontsize_x)
        if fig_save_flag == True:
            plt.savefig('{}/tool_ES_avg{}.pdf'.format(fig_save_dir,figtail), pad_inches = 0.1, bbox_inches = 'tight', dpi=500)
        plt.show()
    
    if verbose:
        for t_i, tmp_tool in enumerate(tool_list):
            print('{}: {}'.format(tmp_tool, all_avg_list.loc[tmp_tool]))

            

def plot_es_workflow(d_rat_df, tool_res_dic, ct_distype_sr, tool_name_trans_dic, fig_dir, pkl_save_dir, 
                     top_per=0.1, origin_gsea=True, drat_pval_cutoff=0.01, pval_cutoff_flag=False,
                     fig_save_flag=True, pkl_save_flag=True, color_diff=True,
                    verbose=True, save_tail_name=''):
    

    ES_dic, cumsum_dic = generate_ES_dic_with_cumsum(d_rat_df, tool_res_dic, ct_distype_sr, origin=origin_gsea, top_per=top_per, pval_cutoff_flag=pval_cutoff_flag, drat_pval_cutoff=drat_pval_cutoff)

    plot_ES_scatter(ES_dic, fig_dir, tool_name_trans_dic, 
                    avg=True,s=100,fontsize=14, 
                    fig_save_flag=fig_save_flag, color_diff=color_diff, verbose=verbose)
    
    if pkl_save_flag == True:
        
        if not os.path.exists(pkl_save_dir):
            os.makedirs(pkl_save_dir)
        
        with open('{}/tool_ES_dic_sinkhorn2{}.pkl'.format(pkl_save_dir, save_tail_name),'wb') as f:
            pkl.dump(ES_dic, f)
        with open('{}/tool_cumsum_dic_sinkhorn2{}.pkl'.format(pkl_save_dir, save_tail_name),'wb') as f:
            pkl.dump(cumsum_dic, f)

            
def plot_simulate_es_workflow(d_rat_df, tool_res_dic, target_ct_pair_ip_dic, ct_distype_sr, 
                              tool_name_trans_dic, fig_dir, pkl_save_dir, 
                              origin_ls = False, verbose=False, origin_gsea=True,
                              top_per=0.1, drat_pval_cutoff=0.01, pval_cutoff_flag=False, color_diff=False, fig_save_flag=True, pkl_save_flag=True, save_tail_name=''):
    

    ES_dic, cumsum_dic = generate_simulate_ES_dic_with_cumsum(d_rat_df, tool_res_dic, 
                                                              target_ct_pair_ip_dic, ct_distype_sr,
                                                             top_per=top_per, origin_ls=origin_ls, origin_gsea=origin_gsea, pval_cutoff_flag=pval_cutoff_flag, drat_pval_cutoff=drat_pval_cutoff)

    plot_ES_scatter(ES_dic, fig_dir, tool_name_trans_dic, 
                    avg=True,s=100,fontsize=14, 
                    fig_save_flag=fig_save_flag, figtail=save_tail_name, verbose=verbose, 
                    color_diff=color_diff)
    
    if pkl_save_flag == True: 
        
        if not os.path.exists(pkl_save_dir):
            os.makedirs(pkl_save_dir)
        
        with open('{}/tool_ES_dic_sinkhorn2{}.pkl'.format(pkl_save_dir, save_tail_name),'wb') as f:
            pkl.dump(ES_dic, f)
        with open('{}/tool_cumsum_dic_sinkhorn2{}.pkl'.format(pkl_save_dir, save_tail_name),'wb') as f:
            pkl.dump(cumsum_dic, f)


            
def plot_ES_multi_scatter(ES_dic_list, tool_name_trans_dic, fig_save_dir,
                    h=6, w=2.5, fontsize=14, fontsize_x=12, 
                          title_flag=True, figtail='', 
                    fig_save_flag=False, verbose=False, color_diff=False):
    
    # plot boxplot for multi-sample
    plot_list = []
    for ES_dic in ES_dic_list:
        tool_list = sorted(list(ES_dic.keys()))
        tool_display_list = [tool_name_trans_dic[tool] for tool in tool_list]

        tmp_all_avg_list = []
        for tool in tool_list:
            if len(ES_dic[tool]['short']) > 0 and len(ES_dic[tool]['long']) > 0:
                tmp_all_avg_list.append((np.average(ES_dic[tool]['short'])+np.average(ES_dic[tool]['long']))/2)
            elif len(ES_dic[tool]['short']) == 0:
                tmp_all_avg_list.append(np.average(ES_dic[tool]['long']))
            elif len(ES_dic[tool]['long']) == 0:
                tmp_all_avg_list.append(np.average(ES_dic[tool]['short']))

        plot_list.append(tmp_all_avg_list)
    
    plot_df = pd.DataFrame(plot_list, columns = tool_display_list).T
    
    avg_es = plot_df.apply(np.nanmedian, axis=1)
    avg_rank = avg_es.rank(ascending=False)
    avg_rank = avg_rank.sort_values(ascending=False)
    rank_list = list(avg_rank.index)
    
    plot_df = plot_df.T.loc[:,rank_list]
        
    plt.figure(figsize=(w,h))
    
    if color_diff:
        color_list = ['red' if round(tmp_es,6) > -1 else 'gray' for tmp_es in list(avg_es.loc[rank_list])]
    else:
        color_list = ['red']*len(rank_list)
        
    bp = plt.boxplot(plot_df, vert=False, patch_artist=True)
    
    for bi, box in enumerate(bp['boxes']):
        box.set(edgecolor='black')
        box.set(facecolor = [1,128/255,0,0.5] )
    for mi, median in enumerate(bp['medians']):
        median.set(color =color_list[mi],linewidth = 1)
    
    plt.yticks([i+1 for i in range(len(rank_list))], rank_list, fontsize=fontsize)
    if title_flag:
        plt.title('DES', fontsize=fontsize+2, pad=10, fontweight='bold')
    plt.xticks(fontsize=fontsize_x)
    
    if fig_save_flag:
        if not os.path.exists(fig_save_dir):
            os.makedirs(fig_save_dir)
        
        plt.savefig('{}/tool_ES_avg_boxplot{}.pdf'.format(fig_save_dir,figtail), pad_inches = 0.1, bbox_inches = 'tight', dpi=500)
    plt.show()
    
    if verbose:
        print(avg_es)
                        
            
    
def split_mulit_ip_d_rat(d_rat_df):
    select_index = ['+' in ip for ip in list(d_rat_df['ip'])]
    tmp_d_rat_df = d_rat_df.loc[select_index,:]
    tmp_d_rat_df_copy = d_rat_df.loc[select_index,:]

    tmp_d_rat_df['ip'] = tmp_d_rat_df['ip'].apply(lambda x: '{} - {}'.format(x.split(' - ')[0], x.split(' - ')[1].split('+')[0][1:]))
    tmp_d_rat_df_copy['ip'] = tmp_d_rat_df_copy['ip'].apply(lambda x: '{} - {}'.format(x.split(' - ')[0], x.split(' - ')[1].split('+')[1][:-1]))
    d_rat_df_mix = d_rat_df.append([tmp_d_rat_df, tmp_d_rat_df_copy])
    d_rat_df_mix.drop_duplicates(inplace = True)
    
    return d_rat_df_mix


def generate_cip_level_dic(tool_res_dic,
                          common_level=None):
    '''
    find common interactions in different levels (common in 2 tools, 3 tools, ...)
    common_level indicate the common tool number
    '''
    
    # find all the common celltypes
    # only cal common interactions in the celltypes existed in both tools' results
    comm_ct_set = ()
    for n, tool in enumerate(tool_res_dic):
        if n == 0:
            comm_ct_set = set(list(tool_res_dic[tool].keys()))
        else:
            comm_ct_set = comm_ct_set.intersection(set(list(tool_res_dic[tool].keys())))

    # find common ips in different levels in each celltype
    # (common in at least 2 tools , 3 tools ,..., et al)
    if common_level == None:
        # no specific level defined
        common_levels = [l for l in range(2,len(list(tool_res_dic.keys())))]
    else:
        common_levels = [common_level]

    cip_level_dic = {}

    for c_level in common_levels:
        tmp_cip_dic = {}
        cip_level_dic.update({c_level:{}})
        for ct in comm_ct_set:
            tmp_ip_list = []
            for tool in tool_res_dic:
                tmp_ips = tool_res_dic[tool][ct]
                # split multi-subunit ips to single
                for ip in tmp_ips:
                    if ip.count('+') > 0:
                        lig, recs = ip.split(' - ')
                        recs = recs[1:-1].split('+')
                        single_ips = ['{} - {}'.format(lig, rec) for rec in recs]
                        tmp_ips += single_ips
                tmp_ips = list(set(tmp_ips))
                tmp_ip_list += tmp_ips
                
            tmp_ip_count = pd.Series(tmp_ip_list).value_counts()
            tmp_ip_count = tmp_ip_count.loc[tmp_ip_count >= c_level]
            tmp_ip_count = tmp_ip_count.sort_values(ascending=False)
            cip_level_dic[c_level].update({
                ct: tmp_ip_count
            })
    
    return cip_level_dic
       

def plot_f1_boxplot(tool_res_dic, cip_level_dic, fig_save_dir, 
                    tool_name_trans_dic=None, 
                    h=6, w=2.5, fontsize=10, color_diff=True,
                    fig_tail='', save_fig=True, return_f1=False):
        
    total_f1_list = []
    tool_list = list(tool_res_dic.keys())
    
    keep_ct_list = [] # record how many celltpe pairs have common interactions
    for ct in cip_level_dic:
        f1_list = []
        cip_list = list(cip_level_dic[ct].index)
        if len(cip_list) == 0:
            # if no common interactions in this celltpe pair, skip
            continue
        else:
            keep_ct_list.append(ct)
        for tool in tool_list:
            tool_res_ip = tool_res_dic[tool][ct]
            for ip in tool_res_ip:
                if ip.count('+') > 0:
                    lig, recs = ip.split(' - ')
                    recs = recs[1:-1].split('+')
                    single_ips = ['{} - {}'.format(lig, rec) for rec in recs]
                    tool_res_ip += single_ips
            tool_res_ip = list(set(tool_res_ip))
    
            ### cal the TP and FN
            if len(tool_res_ip) > 0:
                TP = len(set(tool_res_ip).intersection(set(cip_list)))
                FN = len([ip for ip in cip_list if ip not in tool_res_ip])
                FP = len([ip for ip in tool_res_ip if ip not in cip_list])
                if TP+FP > 0 and TP > 0:
                    precision = TP/(TP+FP)
                    recall = TP/(TP+FN)
                    f1 = 2*precision*recall/(precision+recall)                
                else:
                    f1=0
            else:
                f1 = 0
            f1_list.append(f1)
        total_f1_list.append(f1_list)
        
    f1_df = pd.DataFrame(total_f1_list, index=keep_ct_list, columns=tool_list)
    
    median_f1 = f1_df.apply(np.nanmedian, axis=0).sort_values()
    mean_f1 = f1_df.apply(np.nanmean, axis=0).sort_values()
    # if multi tools have a median of 0
    # rerank these tools by mean
    if median_f1.loc[median_f1 == 0].shape[0] > 0:
        zero_tools = list(median_f1.loc[median_f1 == 0].index)
        replace_rank = [t for t in list(mean_f1.index) if t in zero_tools]
        median_f1.index = replace_rank + list(median_f1.loc[median_f1 != 0].index)
    f1_df = f1_df.loc[:,median_f1.index]
    rank_list = list(median_f1.index)
    
#     f1_dic = {}
#     for tool in rank_list:
#         f1_dic.update({
#             tool: [f1 for f1 in f1_df.loc[:,tool] if not np.isnan(f1)]
#         })
        
    plt.figure(figsize=(w,h))

    bp = plt.boxplot(f1_df, vert=False, patch_artist=True)
    
    if color_diff:
        color_list = ['red' if round(x,2) > 0 else 'gray' for x in median_f1]
    else:
        color_list = ['red']*len(rank_list)

    for bi, box in enumerate(bp['boxes']):
        box.set(edgecolor='black')
        box.set(facecolor = [0,191/255,1,0.5] )
    for mi, median in enumerate(bp['medians']):
        median.set(color =color_list[mi],linewidth = 1)

    if tool_name_trans_dic != None:
        plt.yticks([i+1 for i in range(len(rank_list))], [tool_name_trans_dic[t] for t in rank_list], fontsize=fontsize)
    else:
        plt.yticks([i+1 for i in range(len(rank_list))], rank_list, fontsize=fontsize)
    plt.title('F1 Score', fontsize=fontsize+2, pad=10, fontweight='bold')
    plt.xticks(fontsize=12)
    if save_fig:
        plt.savefig('{}/tool_F1_boxplot{}.pdf'.format(fig_save_dir, fig_tail), pad_inches = 0.1, bbox_inches = 'tight', dpi=500)
    plt.show()
    
    if return_f1:
        return f1_df

    
    
def plot_f1(tool_res_dic, cip_level_dic, c_level=3, tool_name_trans_dic=None, fig_save_dir='./', 
            order_tool_list=None, h=6, w=2.5, s=30, marker='>', color='navy', fontsize=10, fig_tail = '', save_fig = True):
        
    test_cip_level_dic = cip_level_dic[c_level]
    total_f1_list = []
    tool_list = list(tool_res_dic.keys())
    
    keep_ct_list = []
    for ct in test_cip_level_dic:
        f1_list = []
        cip_list = list(test_cip_level_dic[ct].index)
        if len(cip_list) == 0:
            continue
        else:
            keep_ct_list.append(ct)
        for tool in tool_res_dic:
            tool_res_ip = tool_res_dic[tool][ct]

            if tool in ['cc','cpdb','icellnet']:
                ## split multi subunit 
                single_ip_list = [ip for ip in tool_res_ip if ip.count('+') == 0]
                multi_ip_list = [ip for ip in tool_res_ip if ip.count('+') > 0]
                tmp_ip_sub = []
                for n, multi_ip in enumerate(multi_ip_list): 
                    multi_ip = multi_ip.replace('(','').replace(')','').replace('+',' - ').split(' - ')
                    multi_ip_sub = ['{} - {}'.format(multi_ip[0],multi_ip[si]) for si in range(1,len(multi_ip))]
                    for ip_sub in multi_ip_sub:
                        if ip_sub not in tmp_ip_sub:
                            tmp_ip_sub.append(ip_sub)
                single_ip_list += tmp_ip_sub
                tool_res_ip = list(set(single_ip_list))
    
            ### cal the TP and FN
            if len(tool_res_ip) > 0:
                TP = len(set(tool_res_ip).intersection(set(cip_list)))
                FN = len([ip for ip in cip_list if ip not in tool_res_ip])
                FP = len([ip for ip in tool_res_ip if ip not in cip_list])
                if TP+FP > 0 and TP > 0:
                    precision = TP/(TP+FP)
                    recall = TP/(TP+FN)
                    f1 = 2*precision*recall/(precision+recall)                
                else:
                    f1=0
            else:
                f1 = 0
            f1_list.append(f1)
        total_f1_list.append(f1_list)
        
    f1_df = pd.DataFrame(total_f1_list, index=keep_ct_list, columns=tool_list)
    
    plot_list = list(np.average(f1_df,axis=0))
    if tool_name_trans_dic != None:
        plot_list = pd.Series(plot_list, index=[tool_name_trans_dic[tool] for tool in list(f1_df.columns)])
    else:
        plot_list = pd.Series(plot_list, index=list(f1_df.columns))
    if order_tool_list == None:
        plot_list = plot_list.sort_values()
    else:
        plot_list = plot_list.loc[order_tool_list]
        
    plt.figure(figsize=(w,h))
    plt.scatter(plot_list, range(plot_list.shape[0]), color=color, marker=marker, s=s)
    plt.yticks(range(plot_list.shape[0]), list(plot_list.index), fontsize=fontsize)
    plt.title('F1 Score', fontsize=fontsize+2, pad=10, fontweight='bold')
    plt.xticks(fontsize=12)
#     plt.xlim(0.5,1)
    if save_fig == True:
        plt.savefig('{}/tool_F1_avg{}.pdf'.format(fig_save_dir, fig_tail), pad_inches = 0.1, bbox_inches = 'tight', dpi=500)
    plt.show()
    


def generate_f1(tool_res_dic, cip_level_dic, c_level, tool_name_trans_dic):
        
    test_cip_level_dic = cip_level_dic[c_level]
    total_f1_list = []
    tool_list = list(tool_res_dic.keys())
    
    keep_ct_list = []
    for ct in test_cip_level_dic:
        f1_list = []
        cip_list = list(test_cip_level_dic[ct].index)
        if len(cip_list) == 0:
            continue
        else:
            keep_ct_list.append(ct)
        for tool in tool_res_dic:
            tool_res_ip = tool_res_dic[tool][ct]

            if tool in ['cc','cpdb','icellnet']:
                ## split multi subunit 
                single_ip_list = [ip for ip in tool_res_ip if ip.count('+') == 0]
                multi_ip_list = [ip for ip in tool_res_ip if ip.count('+') > 0]
                tmp_ip_sub = []
                for n, multi_ip in enumerate(multi_ip_list): 
                    multi_ip = multi_ip.replace('(','').replace(')','').replace('+',' - ').split(' - ')
                    multi_ip_sub = ['{} - {}'.format(multi_ip[0],multi_ip[si]) for si in range(1,len(multi_ip))]
                    for ip_sub in multi_ip_sub:
                        if ip_sub not in tmp_ip_sub:
                            tmp_ip_sub.append(ip_sub)
                single_ip_list += tmp_ip_sub
                tool_res_ip = list(set(single_ip_list))
    
            ### cal the TP and FN
            if len(tool_res_ip) > 0:
                TP = len(set(tool_res_ip).intersection(set(cip_list)))
                FN = len([ip for ip in cip_list if ip not in tool_res_ip])
                FP = len([ip for ip in tool_res_ip if ip not in cip_list])
                if TP+FP > 0 and TP > 0:
                    precision = TP/(TP+FP)
                    recall = TP/(TP+FN)
                    f1 = 2*precision*recall/(precision+recall)                
                else:
                    f1=0
            else:
                f1 = 0
            f1_list.append(f1)
        total_f1_list.append(f1_list)
        
    f1_df = pd.DataFrame(total_f1_list, index=keep_ct_list, columns=tool_list)
    
    plot_list = list(np.average(f1_df,axis=0))
    plot_list = pd.Series(plot_list, index=[tool_name_trans_dic[tool] for tool in list(f1_df.columns)])
    plot_list = plot_list.sort_values()
    
    return plot_list


def plot_common_metrics_boxplot(tool_res_dic, cip_level_dic, fig_save_dir, 
                    tool_name_trans_dic=None, 
                    h=6, w=2.5, fontsize=10, color_diff=True,
                    fig_tail='', save_fig=True):
        
    total_f1_list, total_precision_list, total_recall_list = [], [], []
    tool_list = list(tool_res_dic.keys())
    
    keep_ct_list = [] # record how many celltpe pairs have common interactions
    for ct in cip_level_dic:
        f1_list, precision_list, recall_list = [], [], []
        cip_list = list(cip_level_dic[ct].index)
        if len(cip_list) == 0:
            # if no common interactions in this celltpe pair, skip
            continue
        else:
            keep_ct_list.append(ct)
        for tool in tool_list:
            tool_res_ip = tool_res_dic[tool][ct]
            for ip in tool_res_ip:
                if ip.count('+') > 0:
                    lig, recs = ip.split(' - ')
                    recs = recs[1:-1].split('+')
                    single_ips = ['{} - {}'.format(lig, rec) for rec in recs]
                    tool_res_ip += single_ips
            tool_res_ip = list(set(tool_res_ip))
    
            ### cal the TP and FN
            if len(tool_res_ip) > 0:
                TP = len(set(tool_res_ip).intersection(set(cip_list)))
                FN = len([ip for ip in cip_list if ip not in tool_res_ip])
                FP = len([ip for ip in tool_res_ip if ip not in cip_list])
                recall = TP/(TP+FN)

                if TP+FP > 0 and TP > 0:
                    precision = TP/(TP+FP)
                    f1 = 2*precision*recall/(precision+recall)                
                else:
                    f1=0
                    precision=0
            else:
                f1 = 0
                precision=0
                recall=0
            f1_list.append(f1)
            precision_list.append(precision)
            recall_list.append(recall)
        total_f1_list.append(f1_list)
        total_precision_list.append(precision_list)
        total_recall_list.append(recall_list)
        
    metric_items = ['F1 Score', 'Precision', 'Recall']
    for n_plot, plot_sources in enumerate([total_f1_list, total_precision_list, total_recall_list]):
        plot_df = pd.DataFrame(plot_sources, index=keep_ct_list, columns=tool_list)

        median_metric = plot_df.apply(np.nanmedian, axis=0).sort_values()
        mean_metric = plot_df.apply(np.nanmean, axis=0).sort_values()
        # if multi tools have a median of 0
        # rerank these tools by mean
        if median_metric.loc[median_metric == 0].shape[0] > 0:
            zero_tools = list(median_metric.loc[median_metric == 0].index)
            replace_rank = [t for t in list(mean_metric.index) if t in zero_tools]
            median_metric.index = replace_rank + list(median_metric.loc[median_metric != 0].index)    
        print(median_metric)
        
        # save the rank in F1 -> rank precision & recall
        if n_plot == 0:
            plot_df = plot_df.loc[:,median_metric.index]
            rank_list = list(median_metric.index)
        else:
            plot_df = plot_df.loc[:,rank_list]
            
        
        ########
        # plot #
        ########
        plt.figure(figsize=(w,h))

        bp = plt.boxplot(plot_df, vert=False, patch_artist=True)

        if color_diff:
            color_list = ['red' if round(x,2) > 0 else 'gray' for x in median_metric]
        else:
            color_list = ['red']*len(rank_list)

        for bi, box in enumerate(bp['boxes']):
            box.set(edgecolor='black')
            box.set(facecolor = [0,191/255,1,0.5] )
        for mi, median in enumerate(bp['medians']):
            median.set(color =color_list[mi],linewidth = 1)

        if tool_name_trans_dic != None:
            plt.yticks([i+1 for i in range(len(rank_list))], [tool_name_trans_dic[t] for t in rank_list], fontsize=fontsize)
        else:
            plt.yticks([i+1 for i in range(len(rank_list))], rank_list, fontsize=fontsize)
        plt.title('{}\nRank'.format(metric_items[n_plot]), fontsize=fontsize+2, pad=10, fontweight='bold')
        plt.xticks(fontsize=12)
        if save_fig:
            plt.savefig('{}/tool_{}_boxplot{}.pdf'.format(fig_save_dir, metric_items[n_plot].replace(' Score',''), fig_tail), pad_inches = 0.1, bbox_inches = 'tight', dpi=500)
        plt.show()
        

def plot_common_metrics_avg_boxplot(tool_res_dic, cip_level_dic, fig_save_dir, 
                    tool_name_trans_dic=None, 
                    h=6, w=2.5, fontsize=10, fontsize_x=12, color_diff=True,
                    fig_tail='', title_flag=True, save_fig=True):
        
    '''
    first avg metrics in each dataset, then boxplot
    '''
    
    total_f1_list, total_precision_list, total_recall_list = [], [], []
    tool_list = list(tool_res_dic.keys())
    
    keep_ct_list = [] # record how many celltpe pairs have common interactions
    for ct in cip_level_dic:
        f1_list, precision_list, recall_list = [], [], []
        cip_list = list(cip_level_dic[ct].index)
        if len(cip_list) == 0:
            # if no common interactions in this celltpe pair, skip
            continue
        else:
            keep_ct_list.append(ct)
        for tool in tool_list:
            tool_res_ip = tool_res_dic[tool][ct]
            for ip in tool_res_ip:
                if ip.count('+') > 0:
                    lig, recs = ip.split(' - ')
                    recs = recs[1:-1].split('+')
                    single_ips = ['{} - {}'.format(lig, rec) for rec in recs]
                    tool_res_ip += single_ips
            tool_res_ip = list(set(tool_res_ip))
    
            ### cal the TP and FN
            if len(tool_res_ip) > 0:
                TP = len(set(tool_res_ip).intersection(set(cip_list)))
                FN = len([ip for ip in cip_list if ip not in tool_res_ip])
                FP = len([ip for ip in tool_res_ip if ip not in cip_list])
                recall = TP/(TP+FN)

                if TP+FP > 0 and TP > 0:
                    precision = TP/(TP+FP)
                    f1 = 2*precision*recall/(precision+recall)                
                else:
                    f1=0
                    precision=0
            else:
                f1 = 0
                precision=0
                recall=0
            f1_list.append(f1)
            precision_list.append(precision)
            recall_list.append(recall)
        total_f1_list.append(f1_list)
        total_precision_list.append(precision_list)
        total_recall_list.append(recall_list)
        
    metric_items = ['F1 Score', 'Precision', 'Recall']
    for n_plot, plot_sources in enumerate([total_f1_list, total_precision_list, total_recall_list]):
        plot_df = pd.DataFrame(plot_sources, index=keep_ct_list, columns=tool_list)
        
        # cal avg by tools
        plot_df['group'] = [bar.split('_')[0] for bar in list(plot_df.index)]
        avg_plot_df = []
        for n_group, sub_df in plot_df.groupby('group'):
            avg_plot_df.append(sub_df.iloc[:,:-1].apply(np.nanmean, axis=0))
        avg_plot_df = pd.DataFrame(avg_plot_df)

        median_metric = avg_plot_df.apply(np.nanmedian, axis=0).sort_values()
        mean_metric = avg_plot_df.apply(np.nanmean, axis=0).sort_values()
        # if multi tools have a median of 0
        # rerank these tools by mean
        if median_metric.loc[median_metric == 0].shape[0] > 0:
            zero_tools = list(median_metric.loc[median_metric == 0].index)
            replace_rank = [t for t in list(mean_metric.index) if t in zero_tools]
            median_metric.index = replace_rank + list(median_metric.loc[median_metric != 0].index)
        print(median_metric)
            
        # save the rank in F1 -> rank precision & recall
        if n_plot == 0:
            avg_plot_df = avg_plot_df.loc[:,median_metric.index]
            rank_list = list(median_metric.index)
        else:
            avg_plot_df = avg_plot_df.loc[:,rank_list]
            
        
        ########
        # plot #
        ########
        plt.figure(figsize=(w,h))

        bp = plt.boxplot(avg_plot_df, vert=False, patch_artist=True)

        if color_diff:
            color_list = ['red' if round(x,2) > 0 else 'gray' for x in median_metric]
        else:
            color_list = ['red']*len(rank_list)

        for bi, box in enumerate(bp['boxes']):
            box.set(edgecolor='black')
            box.set(facecolor = [0,191/255,1,0.5] )
        for mi, median in enumerate(bp['medians']):
            median.set(color =color_list[mi],linewidth = 1)

        if tool_name_trans_dic != None:
            plt.yticks([i+1 for i in range(len(rank_list))], [tool_name_trans_dic[t] for t in rank_list], fontsize=fontsize)
        else:
            plt.yticks([i+1 for i in range(len(rank_list))], rank_list, fontsize=fontsize)
        if title_flag:
            plt.title('{}'.format(metric_items[n_plot]), fontsize=fontsize+2, pad=10, fontweight='bold')
        plt.xticks(fontsize=fontsize_x)
        if save_fig:
            plt.savefig('{}/tool_{}_avg_boxplot{}.pdf'.format(fig_save_dir, metric_items[n_plot].replace(' Score',''), fig_tail), pad_inches = 0.1, bbox_inches = 'tight', dpi=500)
        plt.show()

        
        
        
def plot_common_metrics_scatter(tool_res_dic, cip_level_dic, fig_save_dir='./', 
                    tool_name_trans_dic=None, 
                    h=6, w=2.5, fontsize=10, xticks_font_size=12,
                                title_flag=True,
                    s=30, marker='>', color='navy',
                    fig_tail='', save_fig=True):
        
    total_f1_list, total_precision_list, total_recall_list = [], [], []
    tool_list = list(tool_res_dic.keys())
    
    keep_ct_list = [] # record how many celltpe pairs have common interactions
    for ct in cip_level_dic:
        f1_list, precision_list, recall_list = [], [], []
        cip_list = list(cip_level_dic[ct].index)
        if len(cip_list) == 0:
            # if no common interactions in this celltpe pair, skip
            continue
        else:
            keep_ct_list.append(ct)
        for tool in tool_list:
            tool_res_ip = tool_res_dic[tool][ct]
            for ip in tool_res_ip:
                if ip.count('+') > 0:
                    lig, recs = ip.split(' - ')
                    recs = recs[1:-1].split('+')
                    single_ips = ['{} - {}'.format(lig, rec) for rec in recs]
                    tool_res_ip += single_ips
            tool_res_ip = list(set(tool_res_ip))
    
            ### cal the TP and FN
            if len(tool_res_ip) > 0:
                TP = len(set(tool_res_ip).intersection(set(cip_list)))
                FN = len([ip for ip in cip_list if ip not in tool_res_ip])
                FP = len([ip for ip in tool_res_ip if ip not in cip_list])
                recall = TP/(TP+FN)

                if TP+FP > 0 and TP > 0:
                    precision = TP/(TP+FP)
                    f1 = 2*precision*recall/(precision+recall)                
                else:
                    f1=0
                    precision=0
            else:
                f1 = 0
                precision=0
                recall=0
            f1_list.append(f1)
            precision_list.append(precision)
            recall_list.append(recall)
        total_f1_list.append(f1_list)
        total_precision_list.append(precision_list)
        total_recall_list.append(recall_list)
        
    metric_items = ['F1 Score', 'Precision', 'Recall']
    for n_plot, plot_sources in enumerate([total_f1_list, total_precision_list, total_recall_list]):
        plot_df = pd.DataFrame(plot_sources, index=keep_ct_list, columns=tool_list)
        
        # cal the avg f1, precision, recall across all celltypes

        mean_metric = plot_df.apply(np.nanmean, axis=0).sort_values()
        print(mean_metric)
        
        # save the rank in F1 -> rank precision & recall
        if n_plot == 0:
            rank_list = list(mean_metric.index)
        else:
            mean_metric = mean_metric.loc[rank_list]
            
        
        ########
        # plot #
        ########
        plt.figure(figsize=(w,h))
        plt.scatter(mean_metric, 
                    range(mean_metric.shape[0]), 
                    color=color, 
                    marker=marker, 
                    s=s)
        
        if tool_name_trans_dic != None:
            plt.yticks(range(len(rank_list)), 
                       [tool_name_trans_dic[t] for t in rank_list],
                       fontsize=fontsize)
        else:
            plt.yticks(range(len(rank_list)), 
                       rank_list, fontsize=fontsize)
        if title_flag:
            plt.title('{}\nRank'.format(metric_items[n_plot]), 
                      fontsize=fontsize+2, pad=10, fontweight='bold')
        plt.xticks(fontsize=xticks_font_size)

        if save_fig:
            plt.savefig('{}/tool_{}_scatter{}.pdf'.format(fig_save_dir, metric_items[n_plot].replace(' Score',''), fig_tail), pad_inches = 0.1, bbox_inches = 'tight', dpi=500)
        plt.show()
        



def plot_precision(tool_res_dic, cip_level_dic, c_level, fig_save_dir, tool_name_trans_dic, 
                   order_tool_list=None, h=6, w=2.5, s=30, marker='>', color='navy', fontsize=10, fig_tail='', save_fig = True,):
    
    test_cip_level_dic = cip_level_dic[c_level]
    total_precision_list = []
    tool_list = list(tool_res_dic.keys())
    
    keep_ct_list = []
    for ct in test_cip_level_dic:
        precision_list = []
        cip_list = list(test_cip_level_dic[ct].index)
        if len(cip_list) == 0:
            continue
        else:
            keep_ct_list.append(ct)
        for tool in tool_res_dic:
            tool_res_ip = tool_res_dic[tool][ct]

            if tool in ['cc','cpdb','icellnet']:
                ## split multi subunit 
                single_ip_list = [ip for ip in tool_res_ip if ip.count('+') == 0]
                multi_ip_list = [ip for ip in tool_res_ip if ip.count('+') > 0]
                tmp_ip_sub = []
                for n, multi_ip in enumerate(multi_ip_list): 
                    multi_ip = multi_ip.replace('(','').replace(')','').replace('+',' - ').split(' - ')
                    multi_ip_sub = ['{} - {}'.format(multi_ip[0],multi_ip[si]) for si in range(1,len(multi_ip))]
                    for ip_sub in multi_ip_sub:
                        if ip_sub not in tmp_ip_sub:
                            tmp_ip_sub.append(ip_sub)
                single_ip_list += tmp_ip_sub
                tool_res_ip = list(set(single_ip_list))
    
            ### cal the TP and FN
            if len(tool_res_ip) > 0:
                TP = len(set(tool_res_ip).intersection(set(cip_list)))
                FN = len([ip for ip in cip_list if ip not in tool_res_ip])
                FP = len([ip for ip in tool_res_ip if ip not in cip_list])
                precision = TP/(TP+FP)
            else:
                precision = 0
            precision_list.append(precision)
        total_precision_list.append(precision_list)
        
    precision_df = pd.DataFrame(total_precision_list, index=keep_ct_list, columns=tool_list)
    
    plot_list = list(np.average(precision_df,axis=0))
    plot_list = pd.Series(plot_list, index=[tool_name_trans_dic[tool] for tool in list(precision_df.columns)])
    if order_tool_list == None:
        plot_list = plot_list.sort_values()
    else:
        plot_list = plot_list.loc[order_tool_list]
        
    plt.figure(figsize=(w,h))
    plt.scatter(plot_list, range(plot_list.shape[0]), color=color, marker=marker, s=s)
    plt.yticks(range(plot_list.shape[0]), list(plot_list.index), fontsize=fontsize)
    plt.title('Precision', fontsize=fontsize+2, pad=10, fontweight='bold')
    plt.xticks(fontsize=12)
    if save_fig == True:
        plt.savefig('{}/tool_Precision{}.pdf'.format(fig_save_dir, fig_tail), pad_inches = 0.1, bbox_inches = 'tight', dpi=500)
    plt.show()


def plot_recall(tool_res_dic, cip_level_dic, c_level, fig_save_dir, tool_name_trans_dic, 
                order_tool_list=None, h=6, w=2.5, s=30, marker='>', color='navy', fontsize=10, fig_tail='', save_fig = True):
        
    test_cip_level_dic = cip_level_dic[c_level]
    total_recall_list = []
    tool_list = list(tool_res_dic.keys())
    
    keep_ct_list = []
    for ct in test_cip_level_dic:
        recall_list = []
        cip_list = list(test_cip_level_dic[ct].index)
        if len(cip_list) == 0:
            continue
        else:
            keep_ct_list.append(ct)
        for tool in tool_res_dic:
            tool_res_ip = tool_res_dic[tool][ct]

            if tool in ['cc','cpdb','icellnet']:
                ## split multi subunit 
                single_ip_list = [ip for ip in tool_res_ip if ip.count('+') == 0]
                multi_ip_list = [ip for ip in tool_res_ip if ip.count('+') > 0]
                tmp_ip_sub = []
                for n, multi_ip in enumerate(multi_ip_list): 
                    multi_ip = multi_ip.replace('(','').replace(')','').replace('+',' - ').split(' - ')
                    multi_ip_sub = ['{} - {}'.format(multi_ip[0],multi_ip[si]) for si in range(1,len(multi_ip))]
                    for ip_sub in multi_ip_sub:
                        if ip_sub not in tmp_ip_sub:
                            tmp_ip_sub.append(ip_sub)
                single_ip_list += tmp_ip_sub
                tool_res_ip = list(set(single_ip_list))
    
            ### cal the TP and FN
            TP = len(set(tool_res_ip).intersection(set(cip_list)))
            FN = len([ip for ip in cip_list if ip not in tool_res_ip])
            recall = TP/(TP+FN)
            recall_list.append(recall)
        total_recall_list.append(recall_list)
        
    recall_df = pd.DataFrame(total_recall_list, index=keep_ct_list, columns=tool_list)
    
    plot_list = list(np.average(recall_df,axis=0))
    plot_list = pd.Series(plot_list, index=[tool_name_trans_dic[tool] for tool in list(recall_df.columns)])
    if order_tool_list == None:
        plot_list = plot_list.sort_values()
    else:
        plot_list = plot_list.loc[order_tool_list]
        
    plt.figure(figsize=(w,h))
    plt.scatter(plot_list, range(plot_list.shape[0]), color=color, marker=marker, s=s)
    plt.yticks(range(plot_list.shape[0]), list(plot_list.index), fontsize=fontsize)
    plt.xticks(fontsize=12)
    plt.title('Recall', fontsize=fontsize+2, pad=10, fontweight='bold')
    if save_fig == True:
        plt.savefig('{}/tool_recall{}.pdf'.format(fig_save_dir, fig_tail), pad_inches = 0.1, bbox_inches = 'tight', dpi=500)
    plt.show()

            
def plot_f1_precesion_recall_workflow(tmp_tool_res_dic, tool_name_trans_dic, fig_save_dir, sample, fig_save_flag=True):
    tmp_cip_level_dic = generate_cip_level_dic(tmp_tool_res_dic)

    plot_f1(tmp_tool_res_dic, tmp_cip_level_dic, tool_name_trans_dic = tool_name_trans_dic, c_level=3, fig_save_dir=fig_save_dir, h=6, w=2.5, s=100, marker='>', color='navy', fontsize=14, fig_tail='f1_ranked', save_fig=fig_save_flag)
    tmp_f1_rank = generate_f1(tmp_tool_res_dic, tmp_cip_level_dic, tool_name_trans_dic = tool_name_trans_dic, c_level=3)

    plot_precision(tmp_tool_res_dic, tmp_cip_level_dic, tool_name_trans_dic = tool_name_trans_dic, c_level=3, fig_save_dir=fig_save_dir, h=6, w=2.5, s=100, marker='>', color='navy', fontsize=14, order_tool_list=list(tmp_f1_rank.index), fig_tail='f1_ranked', save_fig=fig_save_flag)
    
    plot_recall(tmp_tool_res_dic, tmp_cip_level_dic, tool_name_trans_dic = tool_name_trans_dic, c_level=3, fig_save_dir=fig_save_dir, h=6, w=2.5, s=100, marker='>', color='navy', fontsize=14, order_tool_list=list(tmp_f1_rank.index), fig_tail='f1_ranked', save_fig=fig_save_flag)

    
    

def extract_tool_res_dic(project_base_dir,
                         tool_list,
                         ct_distype_head_dir='data/pkl',
                         ct_distype_name='ct_distype_sr.pkl',
                         tool_dir_name='tools',
                        cellchat_ip_path='/fs/home/liuzhaoyang/data/cc_ip/cc_ip_all_multi_split_deduplicates.tsv',
                        lr_score=0.6,
                        weight_threshold=0.2,
                        ntop=100,
                        save_dic=True,
                        save_head_dir='evaluation_result/pkl',
                        verbose=True):
    
    
    ################
    # result paths #
    ################
    
    all_cellchat_ips_df = pd.read_csv(cellchat_ip_path, sep='\t', index_col = 0)
    all_cellchat_ips = list(all_cellchat_ips_df.loc[:,'interaction_name_2'])
    
    with open('{}/{}/{}'.format(project_base_dir,ct_distype_head_dir,ct_distype_name), 'rb') as f:
        ct_distype_sr = pkl.load(f)

    ct_list = list(ct_distype_sr.index)
    
    
    #################################
    # extract results for each tool #
    #################################
    
    tool_res_base_dir = '{}/{}'.format(project_base_dir,tool_dir_name)
    tool_res_dic = {}

    for tmp_tool in tool_list:
        res_flag = 1
        
        if tmp_tool == 'cc':
            pval_path_cc = '{}/cc/output/pval_matrix.tsv'.format(tool_res_base_dir)
            name_trans_path_cc = '{}/cc/output/interaction_name_trans.tsv'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_cc(pval_path_cc, ct_list, name_trans_path_cc)
            
        elif tmp_tool == 'cpdb':
            pval_path_cpdb = '{}/cpdb/output/pvalues.txt'.format(tool_res_base_dir)
            name_trans_path_cpdb = '{}/cc/output/interaction_name_trans.tsv'.format(tool_res_base_dir)
            dec_path_cpdb = '{}/cpdb/output/deconvoluted.txt'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_cpdb(pval_path_cpdb, ct_list, name_trans_path_cpdb, dec_path_cpdb)

        elif tmp_tool == 'cpdb_v3':
            pval_path_cpdb = '{}/cpdb_v3/output/pvalues.txt'.format(tool_res_base_dir)
            name_trans_path_cpdb = '{}/cc/output/interaction_name_trans.tsv'.format(tool_res_base_dir)
            dec_path_cpdb = '{}/cpdb_v3/output/deconvoluted.txt'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_cpdb(pval_path_cpdb, ct_list, name_trans_path_cpdb, dec_path_cpdb)

        elif tmp_tool == 'italk':
            res_path_italk = '{}/iTALK/output/LR_result.tsv'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_italk(res_path_italk,ct_list)
            
        elif tmp_tool == 'scr':
            res_path_scr = '{}/SCR/output'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_scr(res_path_scr,ct_list,lr_score=0.6)

        elif tmp_tool == 'cytotalk':
            res_dir_cytotalk = '{}/CytoTalk/output'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_cytotalk(res_dir_cytotalk,ct_list)

        elif tmp_tool == 'natmi':
            res_path_natmi = '{}/NATMI/output/Edges_LRdb_cellchat.csv'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_natmi(res_path_natmi, ct_list, ntop=100, weight_threshold=0.2)

        elif tmp_tool == 'icellnet':
            res_path_icellnet = '{}/icellnet/output'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_icellnet_refine(res_path_icellnet,ct_list)

        elif tmp_tool == 'nichenet':
            res_dir_nichenet = '{}/NicheNet/output'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_nichenet(res_dir_nichenet,ct_list)

        elif tmp_tool == 'scmlnet':
            res_dir_scmlnet = '{}/scmlnet/output'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_scmlnet(res_dir_scmlnet,ct_list)

        elif tmp_tool == 'connectome':
            res_path_connectome = '{}/Connectome/output/weight_sc.tsv'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_connectome(res_path_connectome, ct_list)

        elif tmp_tool == 'cellcall':
            res_path_cellcall = '{}/CellCall/output/lr_exp_log2_scale.tsv'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_cellcall(res_path_cellcall, ct_list)

        elif tmp_tool == 'domino':
            res_dir_domino = '{}/Domino/output/lr_ct_score'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_domino(res_dir_domino, ct_list, all_cellchat_ips)

        elif tmp_tool == 'stlearn':
            res_path_stlearn = '{}/stlearn/output/per_lr_cci_cell_type.tsv'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_stlearn(res_path_stlearn, ct_list)
            
        elif tmp_tool == 'giotto':
            res_path_giotto = '{}/Giotto/output/spatial_scores.tsv'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_giotto(res_path_giotto, ct_list)
            
        elif tmp_tool == 'base_line':
            res_pkl_base_line = '{}/base_line/output/tool_res_dic.pkl'.format(tool_res_base_dir)
            tmp_res_dic = extract_ip_base_line(res_pkl_base_line, ct_list)
            
        else:
            res_flag = 0
            print('unavailable tool name: {}'.format(tmp_tool))
                
        if res_flag == 1:
            tool_res_dic.update({
                tmp_tool: tmp_res_dic
            })
            
            if verbose:
                print('{} done'.format(tmp_tool))
            
    
    ################
    # save results #
    ################
    
    if save_dic:
        project_save_dir = os.path.join(project_base_dir, save_head_dir)
        if not os.path.exists(project_save_dir):
            os.makedirs(project_save_dir)
        with open('{}/tool_res_dic.pkl'.format(project_save_dir), 'wb') as f:
            pkl.dump(tool_res_dic, f)
    
    return tool_res_dic

    