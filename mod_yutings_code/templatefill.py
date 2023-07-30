"""
4/25/2023
Python wrapper script to fill in sbatch job number

Assuming that the template follows this format

#!/bin/bash
#SBATCH -J myjob             # Job name
#SBATCH -o myjob.%j.out      # Name of stdout output file (%j expands to job ID)
#SBATCH -e myjob.%j.err      # Name of stderr output file (%j expands to job ID)
#SBATCH -p normal            # Queue name
#SBATCH -N <node>              # Total number of nodes
#SBATCH -n <mpi>               # Total number of mpi tasks
#SBATCH -t <time>          # Run time (hh:mm:ss)
#SBATCH --dependencies=afterok:<jobid> previousjob
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=<mail>

##TODO: fill in the script path
"""
import sys

def main():
    def slurm_template_fill(id, outscript, previousjob, 
                            inscript='', job='myjob', node='1',
                            mpi='2', time='02:00:00',
                            mail='dhp563@my.utexas.edu'):
        """
        Fill in the template for slurm sbatch script
        
        """
        with open(inscript, 'r') as infile:
            lines = infile.readlines()
        lines[1] = lines[1].replace('myjob', job)
        lines[2] = lines[2].replace('myjob', job)
        lines[3] = lines[3].replace('myjob', job)
        lines[5] = lines[5].replace('<node>', node)
        lines[6] = lines[6].replace('<mpi>', mpi)
        lines[7] = lines[7].replace('<time>', time)
        lines[8] = lines[8].replace('<jobid>', id).replace('previousjob', previousjob)
        lines[10] = lines[10].replace('<mail>', mail)
        with open(outscript, 'w') as outfile:
            outfile.writelines(lines)

    #the inflexible way
    assert len(sys.argv) == 4, "Usage: python templatefill.py <jobid> <script> <jobname> <node> <mpi> <time> <previousjob> <mail>"
    #convert the sys.argv to 
    slurm_template_fill(sys.argv[1], sys.argv[2], sys.argv[3])
    
if __name__ == '__main__':
    main()

        
        
