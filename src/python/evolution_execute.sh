#!/bin/bash                                                                                                                                                                                                        
#Just echoing the numbers as an example
#function_command='echo'

#This is the type of outside call I would make. I just wrote my script (basic_simulation.py) to take an argument at the end
#which it uses as a replicate number to write files. 
#function_command='nohup python basic_simulation.py  ../../Empirical_prot_seqs/Data/structures/1AOE_A.rosetta.pdb'

#So this little loop will run 3000 simulations in batches of 50
#First iterate through the thousands (0, 1, 2)
for i in $(seq 1 2 66); do
    #Then iterate through the hundreds (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
    for j in $(seq 1 20); do
        eval 'nohup python evolution.py paper_data'\$i'.tsv testing.yml' \$j '5000 10 True &' 
    done
    i=$(($i+1))
    for k in $(seq 1 20); do
        eval 'nohup python evolution.py paper_data'\$i'.tsv testing.yml' \$k '5000 10 True &'
    done
    wait #WAIT. This is super important, make sure all those eval calls finish before moving on
        
    #echo "Finished with a batch!"
    #Now iterate through 50-100
    #for k in $(seq 51 100); do  
        #eval 'nohup python evolution.py paper_data'\$i'.tsv testing.yml' \$k '10000 10 True &'
    #done
    #wait #Another wait call to make sure this batch finishes before moving on
        
    #echo "Finished with another batch (which completed a set for one pattern)"
done
