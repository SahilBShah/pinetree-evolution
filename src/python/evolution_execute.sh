#!/bin/bash                                                                                                                                                                                                        

#This loop will run 66 simulations in batches of 25
for i in $(seq 1 2 66); do
    #Iterate through batches of 25
    for j in $(seq 1 25); do
        eval 'nohup python evolution.py paper_data'\$i'.tsv testing.yml' \$j '5000 10 &' 
    done
    i=$(($i+1))
    for k in $(seq 1 25); do
        eval 'nohup python evolution.py paper_data'\$i'.tsv testing.yml' \$k '5000 10 &'
    done
    wait #WAIT. This is super important, make sure all those eval calls finish before moving on
done
