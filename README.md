# Evaluation of cell-cell interaction methods by integrating single-cell RNA sequencing data with spatial information

Recent advances in single-cell RNA-sequencing (scRNA-seq) enable the characterization of cell-cell interactions (CCIs). However, it’s hard to evaluate these methods since no ground truth is provided. Here, we designed a comprehensive workflow to evaluate the accuracy of inferred CCIs from scRNA-seq data based on the expected and observed spatial interaction tendencies. In our workflow, we calculated the metrics of distance enrichment score (DES) and commonly identified interactions to give each tool a comprehensive evaluation.

## CCI evaluation workflow

![image](https://github.com/laputaliu/CCI_evaluation/blob/main/fig/workflow.png)

First, generate known ligand-receptor pairs from CellChatDB, then select short-range and long-range interactions from known pairs for each dataset (top left). Next, perform spatial cell type annotation on ST data coupled with matched scRNA-seq data, and define near and far distributed cell type pairs based on the annotation (medium left). Then feed annotated scRNA-seq data to CCI tools and extract predicted results (bottom left). Finally, evaluate tools’ performances on both distance enrichment score and the metric of commonly identified interactions.

## CCI tools included in this study
- [CellCall](https://github.com/ShellyCoder/cellcall) v0.0.0.9000
- [CellChat](https://github.com/sqjin/CellChat) v1.0.0
- [CellPhoneDB](https://github.com/Teichlab/cellphonedb) v2
- [CellPhoneDB v3](https://github.com/ventolab/CellphoneDB) v3
- [Connectome](https://github.com/msraredon/Connectome) v1.0.1
- [CytoTalk](https://github.com/tanlabcode/CytoTalk) v4.0.11
- [Domino](https://github.com/Chris-Cherry/domino) v0.1.1
- [Giotto](https://github.com/RubD/Giotto) v1.0.4
- [ICELLNET](https://github.com/soumelis-lab/ICELLNET) v0.99.3
- [iTALK](https://github.com/Coolgenome/iTALK) v0.1.0
- [NATMI](https://github.com/forrest-lab/NATMI/)
- [NicheNet](https://github.com/saeyslab/nichenetr) v1.0.0
- [scMLnet](https://github.com/SunXQlab/scMLnet) v0.1.0
- [SingleCellSignalR](https://github.com/SCA-IRCM/SingleCellSignalR_v1) v1.4.0
- [stLearn](https://github.com/BiomedicalMachineLearning/stLearn) v0.4.7

## Datasets included in this study
- <b>Human pancreatic ductal adenocarcinoma (PDAC) dataset</b>
  - The scRNA-seq and ST data of the PDAC dataset can be accessed through GEO under accession number [GSE111672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672).
- <b>Human squamous cell carcinoma (SCC) dataset</b>
  - The scRNA-seq and ST data of the SCC dataset can be accessed through GEO under accession number [GSE144240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144240).
- <b>Mouse cortex dataset</b>
  - The scRNA-seq data of the mouse cortex dataset can be downloaded from [Allen Brain Atlas](https://dx.doi.org/10.1038/nn.4216). The ST data of cortex dataset can be downloaded from [10x Genomics website](https://www.10xgenomics.com/resources/datasets/adult-mouse-brain-ffpe-1-standard-1-3-0).
- <b>Human heart dataset</b>
  - The scRNA-seq and ST data of the human heart dataset can be collected from [Spatial Research website](https://www.spatialresearch.org/resources-published-datasets/doi-10-1016-j-cell-2019-11-025/).
- <b>Human intestinal dataset</b>
  - The scRNA-seq and ST data of the human intestinal dataset can be accessed through GEO under accession numbers [GSE158328](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSE158328) (ST) and [GSE158702](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSE158702) (scRNA-seq)

## Repository

The running scripts of CCI tools and analysis codes used in this study are available in this GitHub repository.
- `prepare_dir_script.py` this Python script is used to prepare the working directories and generate script submitting commands for every CCI tool.
- `/scripts` folder contains all the scripts for running each CCI tool.
- `extract_tool_results.py` this Python script is used to extract all the predicted interactions from each tool's outputs.

<br/>

## Evaluation pipeline
Users can evaluate CCI tools using their own data following the workflow below.

### Step 0: prepare directory and scripts
before we strat, it's a good idea to prepare the directories for storing the results from different tools, we prepared a script for generating all the directories and running scripts for CCI tools inculded in our study. You can simply try it using the example code below, or you can skip this step and prepare directories on your own. 

```
python ./scripts/prepare_dir_script.py --sc_norm ./ST_A3_GSM4797918/data/processed/sc_norm.tsv --sc_count ./ST_A3_GSM4797918/data/processed/sc_counts.tsv --sc_meta ./ST_A3_GSM4797918/data/processed/sc_meta.tsv --deconv ./ST_A3_GSM4797918/data/STRIDE/STRIDE_spot_celltype_frac.txt --st_count ./ST_A3_GSM4797918/data/processed/st_counts.tsv --st_coord ./ST_A3_GSM4797918/data/processed/st_coord.tsv --st_meta ./ST_A3_GSM4797918/data/processed/st_meta.tsv --output_dir ./ST_A3_GSM4797918/tools

```
The running scripts generated here are just for exaples, you can modify them on your own. As for the examples of the generated running scripts, you can refer to `./example_data/ST_A3_GSM4797918/tools/{tool}/script/submit_{tool}.sh` for the script for each tool.

<br/>

### Step 1: compute the d_rat and P-value for LR pairs using ST data
Users need to provie a count file for ST data, a "coordinates" file recording the coordinates of each spot, and a "LR pair" file recording the ligand-receptor interactions. 
- The input "coordinates" file should be a tab-delimited text file with spot barcodes in the first column, the x coordinates of spot in the second column (named as row), and the y coordinated of the spot in the third column (named as col). 
- The "LR pair" file is a LR database in the text format, it should be a tab-delimited text file recording the ligand, receptor and their interactions. You can simply used the "LR pair" file generated from the CellChatDB which is used in our study, you can find it at `./cci_database/cc_ip_multi_split.tsv`. Or you can also use your customized LR database. There are 4 columns needed in the "LR pair", the first column should be the interaction name, you can connect partipated L,R genes with the `_`; next a column named "interaction_name_2" for interaction names in another format (`ligand - receptor` for single subunit, `ligand - (receptor1+receptor2+...)` for multi-subunits); and another two columns recording ligand and receptor gene symbols, for multi-subunit receptors, the subunits need to be seprated by `,`.

<br/>

Example input "coordinates" file:
```
	row	col
CTATCGGGTCTCAACA-1	35	69
TCGAGACCAACACCGT-1	37	55
GATGCGTCCTGCATTC-1	37	57
GCATAGAGCACTCAGG-1	37	59
GCAGATTAGGGATATC-1	37	61
```

<br/>

Example input "LR pair" file (CellChatDB):
```
	pathway_name	annotation	interaction_name_2	ligand	receptor
TGFB1_TGFBR1_TGFBR2	TGFb	Secreted Signaling	TGFB1 - (TGFBR1+TGFBR2)	TGFB1	TGFBR1,TGFBR2
TGFB2_TGFBR1_TGFBR2	TGFb	Secreted Signaling	TGFB2 - (TGFBR1+TGFBR2)	TGFB2	TGFBR1,TGFBR2
TGFB3_TGFBR1_TGFBR2	TGFb	Secreted Signaling	TGFB3 - (TGFBR1+TGFBR2)	TGFB3	TGFBR1,TGFBR2
TGFB1_ACVR1B_TGFBR2	TGFb	Secreted Signaling	TGFB1 - (ACVR1B+TGFBR2)	TGFB1	ACVR1B,TGFBR2
GDF15_TGFBR2	GDF	Secreted Signaling	GDF15 - TGFBR2	GDF15	TGFBR2
GDNF_GFRA1	GDNF	Secreted Signaling	GDNF - GFRA1	GDNF	GFRA1
```

<br/>

when having all these input files prepared, you can submit computing d_rat and P-value of each interaction by using the script we provided, the example command below:
```
nohup python ./scripts/prepare_ip_dis_sinkhorn2.py -c st_count_file -p st_coord_file -o ./ip_dis_sinkhorn2 &
```
Since this step will cost a lot of time, for a quick start, you can use the ip_dis file that we generated before. you can find it in the `./example_data/ST_A3_GSM4797918/data/ip_dis_sinkhorn2/ip_distance_all.tsv`, which computed the d_rat and P-value using the example ST data from sample ST_A3_GSM4797918.


### Step 2: running CCI tools
For this step, becasue of the dependency conflict, we cannot prepare a conda environment for you to run all the tools. So, I am sorry that, you need to prepare the running environment for each tool on your own. But you still can use the scripts generated in the step 0 as an example after having environments prepared. We provided the installation codes and Github links for CCI tools included in our study, you can go to their repositury for the detail.

#### Installation codes
- CellCall
	- for more details, please refer to https://github.com/ShellyCoder/cellcall
```
devtools::install_github("ShellyCoder/cellcall")
```

- CellChat
	- for more details, please refer to https://github.com/sqjin/CellChat
```
devtools::install_github("sqjin/CellChat")
```

- CellPhoneDB & CellPhoneDB v3
	- for more details, please refer to https://github.com/ventolab/CellphoneDB
```
pip install cellphonedb
```

- Connectome
	- for more details, please refer to https://github.com/msraredon/Connectome
```
devtools::install_github('msraredon/Connectome', ref = 'master')
```

- CytoTalk
	- for more details, please refer to https://github.com/tanlabcode/CytoTalk
```
devtools::install_github("tanlabcode/CytoTalk")
```

- Domino
	- for more details, please refer to https://github.com/Chris-Cherry/domino
```
devtools::install_github('Chris-Cherry/domino')
```

- Giotto
	- for more details, please refer to https://github.com/RubD/Giotto
```
remotes::install_github("RubD/Giotto") 
```

- ICELLNET
	- for more details, please refer to https://github.com/soumelis-lab/ICELLNET
```
devtools::install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")
```

- iTALK
	- for more details, please refer to https://github.com/Coolgenome/iTALK
```
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
```

- NATMI
- - for more details, please refer to https://github.com/forrest-lab/NATMI/
```
git clone https://github.com/asrhou/NATMI.git
```

- NicheNet
- - for more details, please refer to https://github.com/saeyslab/nichenetr
```
devtools::install_github("saeyslab/nichenetr")
```

- scMLnet
	- for more details, please refer to https://github.com/SunXQlab/scMLnet
```
devtools::install_github("YUZIXD/scMLnet")
```

- SingleCellSignalR
	- for more details, please refer to https://github.com/SCA-IRCM/SingleCellSignalR_v1
```
devtools::install_github(repo = "https://github.com/SCA-IRCM/SingleCellSignalR_v1", subdir = "SingleCellSignalR")    
```

- stLearn
	- for more details, please refer to https://github.com/BiomedicalMachineLearning/stLearn
```
conda install -c conda-forge stlearn
```

After running these CCI tools, you can use the functions in the `./scripts/extract_tool_results.py` to extract predicted interactions from these CCI tools' output. You can find the example on how to extract tool results from our jupyter notebook file named `extract_results.ipynb` under the `./example_data/ST_A3_GSM4797918/` directory. The original output files of evaluated CCI tools of the sample A3 are also included under the `./example_data/ST_A3_GSM4797918/tools` directory. We stored the extracted tool results in a python dictionary object and dumped it in a pkl file using `pickle` package. You can find the structure of this dictionary in the `including new tools` section. You can find the tool result dictionary of our example dataset at `./example_data/ST_A3_GSM4797918/evaluation_result/pkl/tool_res_dic.pkl` and load it using the code below:
```
import pickle as pkl

with open('./example_data/ST_A3_GSM4797918/evaluation_result/pkl/tool_res_dic.pkl', 'rb') as f:
    tool_res_dic = pkl.load(f)
```

<br/>

#### Including new tools
If users want to including new CCI tools into our evaluation workflow, they only need to provide the results of the addtional tools. The same as what we did above, users can store new tools' results in a python dictionary object and packaged it into a pkl file using pickle. Then, users can provide the path to this result dictionary file to our script and follow the guide in the step 3 to evaluate their new tools. The structure example of the tool result dictionary and the example codes of saving python object using `pickle` are provided below.

Example structure of the tool result dictionary:
- the key is tool name, the value is another dictionary recording predicted interactions in each cell type pair with the key is cell type pair and the value is a list containing the predicted interactions.
```
tool_res_dic = {
    'tool': {
        'celltype_a|celltype_b': ['l1 - r1', 'l2 - r2', 'l3 - (r3_sub1+r3_sub2)'],
        'celltype_a|celltype_c': ['l4 - r4', 'l5 - r5']
    }
}
```

Example codes of saving tool_res_dic using pickle:
```
import pickle as pkl

with open('/your/saving/path/tool_res_dic.pkl', 'wb') as f:
    pkl.dump(tool_res_dic, f)
```

<br/>

### Step 3: calculating distance enrichment score (DES)
After the 2 steps above, then we can finally start to evaluate these CCI tools. You can load the script `./script/extract_tool_reslults.py` for all the functions that we need in this step. And you can find an example on how to evaluate using these functions at `./example_data/ST_A3_GSM4797918/extract_results.ipynb`.

#### Define near/far cell type pairs
Before calculating the DES, we also need to define the near/far cell type pairs based on the average distance between cell types. For doing so, you need to provide an addtional "ST meta" file recording the cell type annotation of the ST data. The spatial cell type annotation can be done by [STRIDE](https://github.com/wanglabtongji/STRIDE). The "ST meta" file should be a tab-delimited text file, the spot barcodes in the first column, the cell types in the second columns.

Example input ST meta file:
```
CTATCGGGTCTCAACA-1	Fibroblasts
TCGAGACCAACACCGT-1	Fibroblasts
GATGCGTCCTGCATTC-1	Fibroblasts
GCATAGAGCACTCAGG-1	Neural
GCAGATTAGGGATATC-1	Neural
```

Having the cell type annotation and the spot coordinates, than we can compute the distance between cell types and define the near/far cell type pairs using the functions `cal_ct_avg_dis()` and `generate_ct_distype()` in the script `extract_tool_results.py`. You can find the details about how to use the example codes below to define near/far cell type pairs in our examples `./example_data/ST_A3_GSM4797918/extract_results.ipynb`

Example codes for defining near/far cell type pairs:
```
# compute the average distance between cell types, using the cell-type annotation and the spot coordinated 
avg_dis_sr = extract_tool_results.cal_ct_avg_dis(meta_df, pos_df)

# define cell type distance type based on the spatial distance 
ct_distype_sr = extract_tool_results.generate_ct_distype(avg_dis_sr)
```

#### Calculating & plotting DES
Through all steps above, we already have all the data that we need for calculating DES: the d_rat file generated in step 1, the tool result dictionary generated in step 2, the cell type pair type generated just now. Then, using the function `plot_es_workflow()` in the script `extract_tool_results.py`, we can compute the DES for each tool in this dataset and plot the DES rank. You can find the details about how to use the example codes below in our examples `./example_data/ST_A3_GSM4797918/extract_results.ipynb`.

Example codes for calculating & plotting DES:
```
import pandas as pd

import sys
sys.path.append('./scripts/')
import extract_tool_results

project_base_dirs = [
    './ST_A3_GSM4797918/'
]

tool_list = [
    'cc', 'cpdb', 'scr', 'natmi', 'icellnet', 'italk',
    'nichenet', 'scmlnet', 'connectome', 'cytotalk', 'cellcall',
    'domino', 'stlearn', 'cpdb_v3', 'giotto', 'base_line'
]

for project_base_dir in project_base_dirs:    
    
    # reading the d_rat file generated in the step 1
    d_rat_df = pd.read_csv('{}/data/ip_dis_sinkhorn2/ip_distance_all.tsv'.format(project_base_dir),
                          sep='\t', index_col = 0)

    # loading the tool_res_dic, and the near/far cell type pair information
    with open('{}/evaluation_result/pkl/tool_res_dic.pkl'.format(project_base_dir), 'rb') as f:
        tool_res_dic = pkl.load(f)
    with open('{}/data/pkl/ct_distype_sr.pkl'.format(project_base_dir), 'rb') as f:
        ct_distype_sr = pkl.load(f)
    
    # calculating & ploting DES
    extract_tool_results.plot_es_workflow(
        d_rat_df, # the d_rat dataframe
        tool_res_dic, # the dictionary recording tool results
        ct_distype_sr, # the pandas Series recording the cell type pair distance type
        '{}/evaluation_result/figure'.format(project_base_dir), # the dir for saving figures
        '{}/evaluation_result/pkl'.format(project_base_dir), # the dir for saving DES scores
        fig_save_flag=True, # whether save the figure
	pkl_save_flag=True, # whether save the DES score
    )    
```

## Data simulation
Still updaing. Coming soon! (^ ^)



