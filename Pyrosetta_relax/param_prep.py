import pyrosetta
pyrosetta.init()

param_topology_file = 'LG.params'

#create
pose = pyrosetta.Pose()
# add paramsblock
rts = pose.conformation().modifiable_residue_type_set_for_conf(pyrosetta.rosetta.core.chemical.FULL_ATOM_t)

topo = open(param_topology_file).read()
buffer = pyrosetta.rosetta.std.stringbuf(topo)
stream = pyrosetta.rosetta.std.istream(buffer)

new = pyrosetta.rosetta.core.chemical.read_topology_file(stream, 'HEM', rts) # no idea what the second argument 
rts.add_base_residue_type(new)
# add new residue
lig = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map( 'HEM' ) )
pose.append_residue_by_jump(lig, 1)
