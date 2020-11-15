#!/bin/bash                                                                                                                                                                                                        

for j in $(seq 1 10); do
    eval 'nohup python evolution.py paper_data13.tsv testing.yml' \$j '5000 10 &'
    eval 'nohup python evolution.py paper_data64.tsv testing.yml' \$j '5000 10 &'
    eval 'nohup python evolution.py paper_data67.tsv testing.yml' \$j '5000 10 &'
    eval 'nohup python evolution.py paper_data79.tsv testing.yml' \$j '5000 10 &'
    eval 'nohup python evolution.py paper_data127.tsv testing.yml' \$j '5000 10 &'
done
wait #WAIT. This is super important, make sure all those eval calls finish before moving on

done
