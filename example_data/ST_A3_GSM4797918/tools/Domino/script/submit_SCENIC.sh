
#!/bin/bash

THREADS=4

WORK_DIR=/fs/home/liuzhaoyang/project/cci_evaluation/human_intestinal/ST_A3_GSM4797918/tools/Domino/SCENIC/scenicdata

cd ${WORK_DIR}

pyscenic grn --num_workers $THREADS -o adjacencies.tsv counts.tsv /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hs_hgnc_curated_tfs.txt


pyscenic ctx adjacencies.tsv /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-500bp-upstream-10species.mc9nr.feather /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-500bp-upstream-7species.mc9nr.feather /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-tss-centered-10kb-10species.mc9nr.feather /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-tss-centered-10kb-7species.mc9nr.feather /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-tss-centered-5kb-10species.mc9nr.feather /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/hg19-tss-centered-5kb-7species.mc9nr.feather --annotations_fname /fs/home/liuzhaoyang/project/cci_evaluation/CCI_tools/Domino/scenicdata_base/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname counts.tsv --mode "dask_multiprocessing" --output regulons.csv --num_workers $THREADS


pyscenic aucell counts.tsv regulons.csv -o auc_mtx.csv --num_workers $THREADS
