#!/bin/bash                                                                                                                                                                                                        

#This loop will run 60 simulations in batches of 50
for i in $(seq 1 2 10); do
    #Iterate through batches of 25
    for j in $(seq 1 6); do
        for k in $seq(1 50); do
            eval 'nohup python evolution.py paper_data'\$i'_arrange'\$j'.tsv testing.yml' \$k '5000 10 &' 
        done    
        wait  
    done
    wait
done
