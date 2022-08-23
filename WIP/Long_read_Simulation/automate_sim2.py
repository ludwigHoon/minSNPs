# STUCK HERE LAST: --args 1 3046 

import time
import os

#id = [4,5,6]
# All ID: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
# Ran: 1
template = """#!/bin/bash
#SBATCH --job-name=ind_{id}_{job_id}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=64:00:00

module load miniconda

source activate r_minsnps

 R CMD BATCH --vanilla "--args {id} {type} {seed}" simulation_run.R
"""
# DONE TILL: 300
for a in range(1,12):
    for type in [1,2]:
        for iii in range(150,300):
            seed = 2050+iii
            new_filename = "ind"+str(iii)+".sh"
            subbed_template = template.replace("{job_id}", str(iii)).replace("{seed}",str(seed)).replace("{id}", str(a)).replace("{type}", str(type))
            #print(seed)
            #print(iii)
            #print(subbed_template)
            with open(new_filename, "w+") as g:
                g.write(subbed_template)
            os.popen('sbatch '+ new_filename).read()

#for file in  independent_sim_*; do echo $file,`sed -n 501p $file` >> line500.txt; done
#for file in  independent_sim_*; do echo $file,`sed -n 196p $file` >> line195.txt; done


#for file in  independent_sim_*_30_10_single*; do echo $file,`tail -n 1 $file` >> sim_summary_30_10_single.txt; done
