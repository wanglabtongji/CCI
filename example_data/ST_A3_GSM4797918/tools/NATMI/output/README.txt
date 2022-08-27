README

ClusterMapping.csv: cell-to-cluster mapping.
Ligands_Receptors_xx.xlsx: information about ligands and receptors in each cell-type/single-cell cluster.
Edges_xx.csv: all ligand-receptor-mediated communications.
LR-pairs_xx folder: all ligand-receptor-mediated communications via a ligand-receptor pair.

Sending cluster: the cluster that expresses the ligand.
Target cluster: the cluster that expresses the receptor.
Ligand/receptor symbol: the official gene symbol of the detected ligand/receptor.
Ligand/receptor-expressing cells: number of cells ligand/receptor is detected in.
Ligand/receptor detection rate: the ratio of cells that expressed the ligand/receptor to total cells in the cluster.
Ligand/receptor average expression value: the average expression level of the ligand/receptor in the cluster.
Ligand/receptor derived specificity of average expression value: the ratio of the average expression level of the ligand/receptor in the cluster to the sum of the average expression levels of the ligand/receptor in every cluster.
Ligand/receptor total expression value: the total expression level of the ligand/receptor in the cluster.
Ligand/receptor derived specificity of total expression value: the ratio of the total expression level of the ligand/receptor in the cluster to the sum of the total expression levels of the ligand/receptor in every cluster.
Edge average expression weight: the product of average expression levels of the ligand and the receptor in the corresponding cluster(s).
Edge average expression derived specificity: the product of the ligand and receptor derived specificity of average expression values.
Edge total expression weight: the product of total expression levels of the ligand and the receptor in the corresponding cluster(s).
Edge total expression derived specificity: the product of the ligand and receptor derived specificity of total expression values.
