#!/bin/bash                                                                                                                                                                                                        

for j in $(seq 1 10); do
    eval 'nohup python evolution.py paper_data1_arrange1.tsv testing.yml' \$j '5000 10 &'
    eval 'nohup python evolution.py paper_data2_arrange2.tsv testing.yml' \$j '5000 10 &'
    eval 'nohup python evolution.py paper_data3_arrange3.tsv testing.yml' \$j '5000 10 &'
    eval 'nohup python evolution.py paper_data4_arrange4.tsv testing.yml' \$j '5000 10 &'
    eval 'nohup python evolution.py paper_data5_arrange5.tsv testing.yml' \$j '5000 10 &'
done
wait #WAIT. This is super important, make sure all those eval calls finish before moving on

done
