#!/bin/bash                                                                                                                                                                                                        

#This will run 10 simulations of a single expression pattern file.
for j in $(seq 1 10); do
    eval 'nohup python evolution.py paper_data1_arrange1.tsv testing.yml' \$j '5000 10 &'
done
wait #WAIT. This is super important, make sure all those eval calls finish before moving on

done
