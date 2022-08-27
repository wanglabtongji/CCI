cat /fs/home/liuzhaoyang/project/cci_evaluation/human_intestinal/ST_A3_GSM4797918/tools/scmlnet/script/scmlnet.input | while read ct_a ct_b count_path meta_path ncore output_dir nohup_path
do
while true
do
sleep 5s
CURRENT=$(ps -ef | grep liuzhao | grep run_scmlnet_multi | wc -l)
if [ ${CURRENT} -lt 8 ];then
screen -x scmlnet -p 0 -X stuff "cd /fs/home/liuzhaoyang/project/cci_evaluation/human_intestinal/ST_A3_GSM4797918 && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_scmlnet_multi.R -a $ct_a -b $ct_b -c $count_path -m $meta_path -n $ncore -o $output_dir > $nohup_path 2>&1 &\n"
break
fi
done
done