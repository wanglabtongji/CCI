import sys
import os
import getopt


usage='Usage:'+sys.argv[0]
usage+='''<Required>[Options]
    <Required>
    -c --sc_norm normalized sc count matrix path
    -r --sc_count sc count matrix path
    -m --sc_meta sc meta path
    -d --deconv deconvolution path, ct fraction matirx
    -s --st_count ST count matrix path
    -p --st_coord ST coord path
    -q --st_meta ST meta path
    -o --output_dir tool output dir
'''

if len(sys.argv) < 2 or not sys.argv[1].startswith('-'):
    sys.exit(usage)

st_meta_path=''

optlist,alist=getopt.getopt(sys.argv[1:],
                            'hc:r:m:d:s:p:q:o:',
                            [
                                'help=',
                                'sc_norm=',
                                'sc_count=',
                                'sc_meta=',
                                'deconv=',
                                'st_count=',
                                'st_coord=',
                                'st_meta=',
                                'output_dir='
                            ])
for opt, arg in optlist:
    if opt in ['-h', '--help']:
        sys.exit(usage)
    elif opt in ['-c', '--sc_norm']:
        count_path = arg
    elif opt in ['-r', '--sc_count']:
        count_raw_path = arg
    elif opt in ['-m', '--sc_meta']:
        meta_path = arg
    elif opt in ['-d', '--deconv']:
        deconv_path = arg
    elif opt in ['-s', '--st_count']:
        st_count_path = arg
    elif opt in ['-p', '--st_coord']:
        st_coord_path = arg
    elif opt in ['-q', '--st_meta']:
        st_meta_path = arg
    elif opt in ['-o', '--output_dir']:
        tool_output_dir = arg


def generate_dir_script(count_path, count_raw_path, meta_path, 
                        deconv_path, st_count_path,st_coord_path, st_meta_path, tool_output_dir):
    '''
    generate the dir and scripts for each tool, (except iTALK)
    - dir : tool - {output, script}
    - script : named as submit_{tool_names}.sh, to run each tool, only need to run ./submit_{tool_names}.sh (maybe)
    
    the screen stuff can not accept the string that is too long, so if the code doesn't work, try to reduce the length of the string in stuff ""
    
    '''
    
    if not os.path.exists(tool_output_dir):
        os.mkdir(tool_output_dir)
    
    if tool_output_dir[-1] == '/':
        tool_output_dir = tool_output_dir[-1]

    ## prepare output dir

    tool_list = ['cc','cpdb','iTALK','SCR','NATMI','icellnet','NicheNet','CytoTalk','scmlnet','Connectome','Domino','CellCall','stlearn','Giotto','cpdb_v3']

    for tool in tool_list:
        try:
            os.mkdir(os.path.join(tool_output_dir, tool))
            os.mkdir(os.path.join(tool_output_dir, tool, 'script'))
            os.mkdir(os.path.join(tool_output_dir, tool, 'output'))
        except FileExistsError:
            continue

    ## modify input path, to avoid the string which send to screen is too long
    tmp_input_list = [count_path, count_raw_path, meta_path, deconv_path, st_count_path, st_coord_path, tool_output_dir]
    data_base_dir = os.path.commonpath(tmp_input_list)
    
    count_path_short = count_path.replace(data_base_dir,'.')
    count_raw_path_short = count_raw_path.replace(data_base_dir,'.')
    meta_path_short = meta_path.replace(data_base_dir,'.')
    deconv_path_short = deconv_path.replace(data_base_dir,'.')
    st_count_path_short = st_count_path.replace(data_base_dir,'.')
    st_coord_path_short = st_coord_path.replace(data_base_dir,'.')
    st_meta_path_short = st_meta_path.replace(data_base_dir,'.')
    tool_output_dir_short = tool_output_dir.replace(data_base_dir,'.')
    
    
    
    
    #########################################
    # generate running script for each tool #
    #########################################
    
    
    
    ############
    # cellchat #
    ############
    

    cc_command = 'screen -x cellchat -p 0 -X stuff "cd {} && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_cc.R -c {} -m {} -o {}/cc/output/ > {}/cc/script/cc_nohup.out 2>&1 &\\n"'.format(data_base_dir, count_path_short, meta_path_short, tool_output_dir_short, tool_output_dir_short)

    with open('{}/cc/script/submit_cc.sh'.format(tool_output_dir),'w') as f:
        f.write(cc_command)
        
    os.system('chmod 777 {}/cc/script/submit_cc.sh'.format(tool_output_dir))

    
    
    
    ###############
    # cellPhoneDB #
    ###############
    
    
    cpdb_command = 'screen -x cpdb -p 0 -X stuff "cd {} && nohup time -v cellphonedb method statistical_analysis {} {} --counts-data gene_name --output-path {}/cpdb/output --threads 4 --database /fs/home/liuzhaoyang/data/cc_ip/cellphonedb_user_2021-01-12-18_19.db > {}/cpdb/script/cpdb_nohup.out 2>&1 &\\n"'.format(data_base_dir, meta_path_short, count_path_short, tool_output_dir_short, tool_output_dir_short)

    with open('{}/cpdb/script/submit_cpdb.sh'.format(tool_output_dir),'w') as f:
        f.write(cpdb_command)
        
    os.system('chmod 777 {}/cpdb/script/submit_cpdb.sh'.format(tool_output_dir))


    
    #########
    # iTALK #
    #########
    
    
    italk_command = 'screen -x italk -p 0 -X stuff "cd {} && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_italk.R -c {} -m {} -o {}/iTALK/output/ > {}/iTALK/script/italk_nohup.out 2>&1 &\\n"'.format(data_base_dir, count_path_short, meta_path_short, tool_output_dir_short, tool_output_dir_short)

    with open('{}/iTALK/script/submit_italk.sh'.format(tool_output_dir),'w') as f:
        f.write(italk_command)
        
    os.system('chmod 777 {}/iTALK/script/submit_italk.sh'.format(tool_output_dir))

    
    
    
    #####################
    # SingleCellSignalR #
    #####################
    
    
    scr_command = 'screen -x scr -p 0 -X stuff "cd {} && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_scr_cutoff_refine.R -c {} -m {} -o {}/SCR/output/ -s 0.6 > {}/SCR/script/scr_nohup.out 2>&1 &\\n"'.format(data_base_dir, count_path_short, meta_path_short, tool_output_dir_short, tool_output_dir_short)

    with open('{}/SCR/script/submit_scr.sh'.format(tool_output_dir),'w') as f:
        f.write(scr_command)
        
    os.system('chmod 777 {}/SCR/script/submit_scr.sh'.format(tool_output_dir))


    
    #########
    # NATMI #
    #########
    

    natmi_command1 = 'sed \'1 iCell\\tAnnotation\' {} > {}/NATMI/meta_natmi.tsv && cd ~/biosoft/NATMI && time -v python ExtractEdges.py --interDB LRdb_cellchat --interSpecies human --emFile {} --annFile {}/NATMI/meta_natmi.tsv --species human --idType symbol --out {}/NATMI/output/'.format(meta_path, tool_output_dir, count_path, tool_output_dir, tool_output_dir)

    with open('{}/NATMI/script/run_natmi.sh'.format(tool_output_dir),'w') as f:
        f.write(natmi_command1)

    natmi_command = 'chmod 777 {}/NATMI/script/run_natmi.sh && screen -x natmi -p 0 -X stuff "cd {}/NATMI/script/ && nohup ./run_natmi.sh > {}/NATMI/script/natmi_nohup.out 2>&1 &\\n"'.format(tool_output_dir, tool_output_dir, tool_output_dir)

    with open('{}/NATMI/script/submit_natmi.sh'.format(tool_output_dir),'w') as f:
        f.write(natmi_command)
        
    os.system('chmod 777 {}/NATMI/script/submit_natmi.sh'.format(tool_output_dir))


    

    ############
    # icellnet #
    ############

    
    icellnet_command = 'screen -x icellnet -p 0 -X stuff "cd {} && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_icellnet_refine.R -c {} -m {} -o {}/icellnet/output/ > {}/icellnet/script/icellnet_nohup.out 2>&1 &\\n"'.format(data_base_dir, count_path_short, meta_path_short, tool_output_dir_short, tool_output_dir_short)

    with open('{}/icellnet/script/submit_icellnet.sh'.format(tool_output_dir),'w') as f:
        f.write(icellnet_command)
        
    os.system('chmod 777 {}/icellnet/script/submit_icellnet.sh'.format(tool_output_dir))

    
    
    ############
    # NicheNet #
    ############
    
    
    nichenet_command = 'screen -x nichenet -p 0 -X stuff "cd {} && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_nichenet.R -c {} -m {} -o {}/NicheNet/output/ > {}/NicheNet/script/nichenet_nohup.out 2>&1 &\\n"'.format(data_base_dir, count_path_short, meta_path_short, tool_output_dir_short, tool_output_dir_short)

    with open('{}/NicheNet/script/submit_nichenet.sh'.format(tool_output_dir),'w') as f:
        f.write(nichenet_command)
        
    os.system('chmod 777 {}/NicheNet/script/submit_nichenet.sh'.format(tool_output_dir))

    
    
    ############
    # CytoTalk #
    ############

    
    ## the prepare_input_cytotalk.py will generate all the files of the input dir and also for the run_cytotalk.sh and cytotalk.input
    cytotalk_command = 'python /fs/home/liuzhaoyang/project/cci_evaluation/scripts/prepare_input_cytotalk.py -c {} -m {} -o {} && cd {}/CytoTalk/script/ && nohup time -v ./run_cytotalk.sh > {}/CytoTalk/script/cytotalk_nohup.out 2>&1 &'.format(count_path, meta_path, tool_output_dir, tool_output_dir, tool_output_dir)

    ## for the real run, using the submit_cytotalk.sh
    with open('{}/CytoTalk/script/submit_cytotalk.sh'.format(tool_output_dir),'w') as f:
        f.write(cytotalk_command)
        
    os.system('chmod 777 {}/CytoTalk/script/submit_cytotalk.sh'.format(tool_output_dir))
    
    ## if it occurs the error when generate the output, then we can use the script/regenerate_cytotalk_result.R


    ###########
    # scmlnet #
    ###########
    
    
    # scmlnet_command = 'screen -x scmlnet -p 0 -X stuff "nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_scmlnet.R -c {} -m {} -o {}/scmlnet/output/ -n {} > {}/scmlnet/script/scmlnet_nohup.out 2>&1 &\\n"'.format(count_raw_path, meta_path, tool_output_dir, 4, tool_output_dir)

    # with open('{}/scmlnet/script/submit_scmlnet.sh'.format(tool_output_dir),'w') as f:
    #     f.write(scmlnet_command)

    # os.system('chmod 777 {}/scmlnet/script/submit_scmlnet.sh'.format(tool_output_dir))

    scmlnet_command = 'sed \'1 iBarcode\tCluster\' {} > {}/scmlnet/scmlnet_meta.tsv && python /fs/home/liuzhaoyang/project/cci_evaluation/scripts/prepare_input_scmlnet.py -c {} -m {}/scmlnet/scmlnet_meta.tsv -o {} -t 1 && cd {}/scmlnet/script/ && nohup time -v ./run_scmlnet.sh > {}/scmlnet/script/scmlnet_nohup.out 2>&1 &'.format(meta_path, tool_output_dir, count_raw_path, tool_output_dir, tool_output_dir, tool_output_dir, tool_output_dir)

    with open('{}/scmlnet/script/submit_scmlnet.sh'.format(tool_output_dir),'w') as f:
        f.write(scmlnet_command)

    os.system('chmod 777 {}/scmlnet/script/submit_scmlnet.sh'.format(tool_output_dir))



    ##############
    # Connectome #
    ##############
    
    connectome_command = 'screen -x connectome -p 0 -X stuff "cd {} && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_connectome.R -c {} -m {} -o {}/Connectome/output/ > {}/Connectome/script/connectome_nohup.out 2>&1 &\\n"'.format(data_base_dir, count_path_short, meta_path_short, tool_output_dir_short, tool_output_dir_short)

    with open('{}/Connectome/script/submit_connectome.sh'.format(tool_output_dir),'w') as f:
        f.write(connectome_command)
        
    os.system('chmod 777 {}/Connectome/script/submit_connectome.sh'.format(tool_output_dir))


    
    ##########
    # Domino #
    ##########
    
    
    # step 1: running SCENIC
    
    if not os.path.exists('{}/Domino/SCENIC'.format(tool_output_dir)):
        os.mkdir(os.path.join(tool_output_dir,'Domino','SCENIC'))
        os.mkdir(os.path.join(tool_output_dir,'Domino','SCENIC','scenicdata'))
        
    # prepare input file for SCENIC
    prepare_scenic_command = "python /fs/home/liuzhaoyang/project/cci_evaluation/scripts/prepare_scenic_input.py -c {} -o {}".format(count_raw_path, os.path.join(tool_output_dir,'Domino','SCENIC','scenicdata','counts.tsv'))
    
    with open('{}/Domino/script/prepare_scenic_input.sh'.format(tool_output_dir),'w') as f:
        f.write(prepare_scenic_command)
        
    os.system('chmod 777 {}/Domino/script/prepare_scenic_input.sh'.format(tool_output_dir))

        
    # SCENIC running script
    scenic_command = '''
#!/bin/bash

THREADS=4

WORK_DIR={}

cd ${{WORK_DIR}}

pyscenic grn --num_workers $THREADS -o adjacencies.tsv counts.tsv /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hs_hgnc_curated_tfs.txt


pyscenic ctx adjacencies.tsv /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-500bp-upstream-10species.mc9nr.feather /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-500bp-upstream-7species.mc9nr.feather /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-tss-centered-10kb-10species.mc9nr.feather /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-tss-centered-10kb-7species.mc9nr.feather /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-tss-centered-5kb-10species.mc9nr.feather /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-tss-centered-5kb-7species.mc9nr.feather --annotations_fname /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname counts.tsv --mode "dask_multiprocessing" --output regulons.csv --num_workers $THREADS


pyscenic aucell counts.tsv regulons.csv -o auc_mtx.csv --num_workers $THREADS
'''.format(os.path.join(tool_output_dir,'Domino','SCENIC','scenicdata'))
    
    with open('{}/Domino/script/submit_SCENIC.sh'.format(tool_output_dir),'w') as f:
        f.write(scenic_command)
        
    os.system('chmod 777 {}/Domino/script/submit_SCENIC.sh'.format(tool_output_dir))
    
    
    # prepare input files for Domino
    domino_prepare_script = "Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/prepare_domino_input.R -c {} -m {} -o {}".format(count_raw_path, meta_path, os.path.join(data_base_dir,'data','processed','sc_res.rds'))
    
    with open('{}/Domino/script/prepare_domino_input.sh'.format(tool_output_dir),'w') as f:
        f.write(domino_prepare_script)
        
    os.system('chmod 777 {}/Domino/script/prepare_domino_input.sh'.format(tool_output_dir))
    
    
    # Domino running script
    domino_script = 'screen -x domino -p 0 -X stuff "cd {} && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_domino.R -r ./data/processed/sc_res.rds -s ./tools/Domino/SCENIC/scenicdata -f n -o ./tools/Domino/output/ > ./tools/Domino/script/Domino_nohup.out 2>&1 &\n"'.format(data_base_dir)
    
    with open('{}/Domino/script/submit_Domino.sh'.format(tool_output_dir),'w') as f:
        f.write(domino_script)
        
    os.system('chmod 777 {}/Domino/script/submit_Domino.sh'.format(tool_output_dir))
    
    
    
    ############
    # CellCall #
    ############
    
    
    cellcall_script = 'screen -x cellcall -p 0 -X stuff "cd {} && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_cellcall.R -c {} -m {} -o {}/CellCall/output > ./tools/CellCall/script/CellCall_nohup.out 2>&1 &\n"'.format(data_base_dir,count_raw_path_short,meta_path_short,tool_output_dir_short)
    
    with open('{}/CellCall/script/submit_cellcall.sh'.format(tool_output_dir),'w') as f:
        f.write(cellcall_script)
        
    os.system('chmod 777 {}/CellCall/script/submit_cellcall.sh'.format(tool_output_dir))
    
    
    
    ###########
    # stLearn #
    ###########
    
    
    stlearn_script = 'screen -x stlearn -p 0 -X stuff "cd {} && nohup time -v python /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_stlearn.py -c {} -p {} -d {} -o ./tools/stlearn/output > ./tools/stlearn/script/stlearn_nohup.out 2>&1 &\n"'.format(data_base_dir, st_count_path_short, st_coord_path_short, deconv_path_short)
    
    with open('{}/stlearn/script/submit_stlearn.sh'.format(tool_output_dir),'w') as f:
        f.write(stlearn_script)
        
    os.system('chmod 777 {}/stlearn/script/submit_stlearn.sh'.format(tool_output_dir))
    
    
    
    ##########
    # Giotto #
    ##########
    
    
    giotto_script = 'screen -x giotto -p 0 -X stuff "cd {} && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_giotto.R -c {} -m {} -p {} -o {}/Giotto/output > {}/Giotto/script/Giotto_nohup.out 2>&1 &\n"'.format(data_base_dir,st_count_path_short,st_meta_path_short, st_coord_path_short, tool_output_dir_short, tool_output_dir_short)
    
    with open('{}/Giotto/script/submit_giotto.sh'.format(tool_output_dir),'w') as f:
        f.write(giotto_script)
        
    os.system('chmod 777 {}/Giotto/script/submit_giotto.sh'.format(tool_output_dir)) 
    
    
    

    
    #########################
    # script for cal ip_dis #
    #########################
    
    if not os.path.exists('{}/data/ip_dis_sinkhorn2'.format(data_base_dir)):
        os.mkdir(os.path.join(data_base_dir,'data','ip_dis_sinkhorn2'))
        
    
    ip_dis_script = 'cd {} && nohup python /fs/home/liuzhaoyang/project/cci_evaluation/scripts/prepare_ip_dis_sinkhorn2.py -c {} -p {} -o {} &'.format(os.path.join(data_base_dir,'data','ip_dis_sinkhorn2'),st_count_path, st_coord_path, os.path.join(data_base_dir,'data','ip_dis_sinkhorn2'))
    
    with open('{}/data/ip_dis_sinkhorn2/submit_ip_dis_sinkhorn2.sh'.format(data_base_dir),'w') as f:
        f.write(ip_dis_script)
        
    os.system('chmod 777 {}/data/ip_dis_sinkhorn2/submit_ip_dis_sinkhorn2.sh'.format(data_base_dir))
    
    
    
    
    


if __name__ == '__main__':
    
#     print('{}\n{}\n{}\n{}\n'.format(count_path, count_raw_path, meta_path, tool_output_dir))
       
    generate_dir_script(count_path, count_raw_path, meta_path, 
                        deconv_path, st_count_path,st_coord_path, st_meta_path, tool_output_dir)


