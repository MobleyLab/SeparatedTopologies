from SeparatedTopologies import make_septop as ms
from SeparatedTopologies import boresch_restraints as br
from SeparatedTopologies import combine_coordinates as ac
from SeparatedTopologies import ligand_files as lf
import os
import parmed as pmd

##########################
###### input needed  #####
##########################

#path to data folder
path = ''
#ligand folders
ligand_A = ''
ligand_B = ''
#directories
compound_A = '%s/%s'%(path, ligand_A)
compound_B = '%s/%s'%(path, ligand_B)
edge_A_B = '%s/edge_%s_%s'%(path, ligand_A, ligand_B)
if not os.path.isdir(edge_A_B):os.mkdir(edge_A_B)

#Three letter code for ligands
lig = 'LIG'

##########################
### Combine .gro files ###
##########################
pdb_A = '%s/complex.pdb'%(compound_A)
pdb_B = '%s/complex.pdb'%(compound_B)
fit_B = '%s/complex_fit.pdb'%compound_B

gro_B = '%s/complex_fit.gro'%compound_B
gro_A = '%s/complex.gro'%(compound_A)
complex = '%s/complex.gro'%edge_A_B

# Align complexed to be able to insert ligand B into complex of ligand A
ac.align_complexes(pdb_A, pdb_B, fit_B)
# Convert .pdb to .gro
ac.pdb2gro(fit_B, gro_B)

#Insert ligand B  coordinates into coordinate file or complex A
complex = ac.combine_ligands_gro(gro_A, gro_B, complex, ligand_A=lig, ligand_B=lig)
# Edit indices
complex = ac.edit_indices(complex,complex)

########################
# Get mol2 ligand file #
########################
#Only if not available yet
mol2_A = '%s/ligand.mol2'%compound_A
mol2_B = '%s/ligand.mol2'%compound_B
lf.ligand_sdf(pdb_A, ligand_A)
lf.ligand_sdf(pdb_B, ligand_B)

##################################
### Compute Boresch restraints ###
##################################

#trajectory or single frame, if a trajectory is provided, that information can help identify stable atoms for Boresch restraints
traj_A = '%s/complex.h5'%(compound_A)
traj_B = '%s/complex.h5'%(compound_B)
# #index of protein atoms has to be zero based
restrain_A, restrain_B, dG_A_off, dG_B_on = br.restrain_ligands(traj_A,traj_B,mol2_A,mol2_B, '%s/boresch_restraints_A_on.itp'%edge_A_B, '%s/boresch_restraints_B_off.itp'%edge_A_B, '%s/boresch_restraints_A.itp'%edge_A_B,'%s/boresch_restraints_B.itp'%edge_A_B,
                                                 ligand_atoms=None, protein_atoms=None, substructure=None,ligand=lig, top_A=gro_A, top_B=gro_B)

###############################
### Make separated topology ###
###############################

#Add ligand B to .top of complex_A

top_A = '%s/complex.top'%(compound_A)
top_B = '%s/complex.top'%(compound_B)
top = '%s/complex.top'%edge_A_B

#Combine topology files
ms.combine_ligands_top(top_A, top_B, top, ligand=lig)

#Name top files
step_1 = '%s/complex_1_scale.top'%edge_A_B
step_2 = '%s/complex_2_scale.top'%edge_A_B

# Generate separated topology files different legs cycle

ms.create_top(top, step_1, 0.5, 'vdwq_scaled-vdw', 'dummy_scaled-vdwq', top_A, top_B, ligand=lig)
br.include_itp_in_top(step_1, 'boresch_restraints_B.itp')
br.include_itp_in_top(step_1, 'boresch_restraints_A_on.itp')

ms.create_top(top, step_2, 0.5, 'scaled-vdw_dummy', 'scaled-vdwq_vdwq', top_A, top_B, ligand=lig)
br.include_itp_in_top(step_2, 'boresch_restraints_B_off.itp')
br.include_itp_in_top(step_2, 'boresch_restraints_A.itp')
