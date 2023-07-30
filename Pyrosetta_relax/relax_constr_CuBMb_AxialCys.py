"""
PyRosetta myoglobin relaxation protocol
"""
from pyrosetta import *
from pyrosetta.rosetta.protocols.simple_moves import SetupMetalsMover
from pyrosetta.teaching import get_fa_scorefxn
from pyrosetta.rosetta.protocols.constraint_movers import ConstraintSetMover 
from rosetta.protocols.relax import FastRelax

##### PARAMS #######
pdbfile='redo3_CuBMb_Cys_axial.pdb'
constraint='constraints.cst'

#Start

init()
#Metal detection

metaldetector = SetupMetalsMover()

# #Fold Tree
# ft = FoldTree()
# ft.add_edge(1,30,-1)
# ft.add_edge(5,31,1)

#Scorefunction with metal binding constr.
scorefxn = get_fa_scorefxn()
scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.metalbinding_constraint, 1.0)

#Constraints
constraints = ConstraintSetMover()
constraints.constraint_file(constraint)
constraints.add_constraints(True)

#Relax function

relax = FastRelax()
relax.set_scorefxn(scorefxn)
relax.constrain_relax_to_start_coords(True)
relax.constrain_coords(True)

#Generation of ensemble
E_relax=[]
for i in range(1,11):
	
	#Create the pose from PDB-File
	globin_Pose = pose_from_pdb(pdbfile)       
	
	#Using the SetupMetalMover
	metaldetector.apply(globin_Pose)        

	# #Set the correct fold tree
	# ZnF_Pose.fold_tree(ft)                                   
	
	#Set constraints
	constraints.apply(globin_Pose)                      
	
	#Relax
	relax.apply(globin_Pose)
	
	#Calculate score of relaxed pose and store it the list E_relax
	E_relax.append(scorefxn(globin_Pose))
	
	#Create PDB-File from relaxed pose
	globin_Pose.dump_pdb('Globin_Pose_relax_{}.pdb'.format(i))
	


#Mean Score, standard deviation and standard error

mean_E_rel = sum(E_relax)/len(E_relax)

A = []
for i in range(0,10):
	a = (E_relax[i]-mean_E_rel)**2
	A.append(a)

std_E_rel = (sum(A) / (len(E_relax)-1))**(0.5)
SEM_E_rel = std_E_rel/((len(E_relax))**(0.5))



#Score of the starting structure

globin_Pose_start = pose_from_pdb(pdbfile)
metaldetector.apply(globin_Pose_start)
constraints.apply(globin_Pose_start)
start_E = scorefxn(globin_Pose_start)
globin_Pose_start.dump_pdb(f'{pdbfile}_pose_start.pdb')


#Outputfile with score of starting structure, mean score of ensemble and score of each relaxed structure 

txtfile = f'{pdbfile}_FastRelax_score.txt'

myfile = open(txtfile, 'w')

myfile.write('start_E = {}\nmean_E_rel = {} ; Std = {} ; SEM = {}\n\nE_rel\n'.format(start_E, mean_E_rel, std_E_rel, SEM_E_rel))

for i in range(0,len(E_relax)):
	myfile.write('{}\n'.format(E_relax[i]))
	
myfile.close()