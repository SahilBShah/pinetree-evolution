#!/bin/bash                                                                                                                                                                                                        

#So this little loop will run 66 simulations in batches of 25
for i in $(seq 1 2 66); do
    #Then iterate through the hundreds (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
    for j in $(seq 26 50); do
        eval 'nohup python evolution.py paper_data'\$i'.tsv testing.yml' \$j '5000 10 True &' 
    done
    i=$(($i+1))
    for k in $(seq 26 50); do
        eval 'nohup python evolution.py paper_data'\$i'.tsv testing.yml' \$k '5000 10 True &'
    done
    wait #WAIT. This is super important, make sure all those eval calls finish before moving on
done
