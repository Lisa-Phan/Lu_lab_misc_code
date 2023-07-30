module load gromacs
gmx=$TACC_GROMACS_BIN/mpi_gmx

sample="1PGB"

$gmx pdb2gmx -f $sample.pdb 
             -o $sample_conf.gro 
             -p $sample_topol.top 
             -ff amber99 
             -water tip3p

$gmx editconf -f $sample_conf.gro 
              -o $sample_solvated.gro
              -bt cubic 
              -d 0.7 

$gmx solvate -cp $sample_solvated.gro
             -cs tip3p
             -p $sample_topol.top
             -o $sample_solvated.gro

$gmx grompp -f sample.mdp 
            -c $sample_solvated.gro
            -p $sample_topol.top
            -o $sample_emin.tpr

$gmx mdrun -v -deffnm $sample_emin


$gmx solvate -cp box.pdb -cs spc216 -o water.pdb -p topol.top
$gmx grompp - 
