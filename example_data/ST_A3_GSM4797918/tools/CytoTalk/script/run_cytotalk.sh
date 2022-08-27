cat /fs/home/liuzhaoyang/project/cci_evaluation/human_intestinal/ST_A3_GSM4797918/tools/CytoTalk/script/cytotalk.input | while read ct_a ct_b input_dir output_dir nohup_path
do
while true
do
CURRENT=$(ps -ef | grep liu | grep -E 'time -v Rscript' | grep -E 'cytotalk' | wc -l)
if [ ${CURRENT} -lt 3 ];then
screen -x cytotalk -p 0 -X stuff "cd /fs/home/liuzhaoyang/project/cci_evaluation/human_intestinal/ST_A3_GSM4797918 && nohup time -v Rscript /fs/home/liuzhaoyang/project/cci_evaluation/scripts/run_cytotalk.R -a $ct_a -b $ct_b -i $input_dir -o $output_dir > $nohup_path 2>&1 &\n"
sleep 20s
break
fi
done
done