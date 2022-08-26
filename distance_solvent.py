import numpy as np
import mdtraj as md
import itertools
from simtk import unit
from scipy.spatial import distance
from openeye import oechem
import parmed as pmd

def distance_restraint(gro,ligand):
    traj = md.load(gro)
    atoms_restrain = []
    for r in traj.topology.residues:
        if r.name == ligand:
            atoms = []
            for a in r.atoms:
                atoms.append(a.index)
            ligand_traj = traj.atom_slice(atoms, inplace=False)
            ligand_coord = ligand_traj.xyz[:,:,:]
            com_ligand = md.compute_center_of_mass(ligand_traj)
            com_distance = []
            #Calculate the distance of ligand atoms to the COM
            for c in ligand_coord[0]:
                d = distance.euclidean(com_ligand[0], c)
                com_distance.append(d)
            # store index with sorted value
            sorted_distance = sorted((x, y) for y, x in enumerate(com_distance))
            #Find ligand heavy atom closest to COM
            index_L1 = sorted_distance[0][1]
            L1 = atoms[index_L1]
            atoms_restrain.append(L1+1)
    unit_cell = traj.unitcell_lengths[-1]
    return atoms_restrain, min(unit_cell)/2
  
def write_itp_restraints(atoms_restrain, dist, forceconst_A, forceconst_B, file):
    """Add dihedral restraints
    Parameters
    ----------
    dih: list
        nested list of ligand atoms selected for restraints
    values: list
        List of values for dihedral
    forceconst_A/B: int
        forceconstant for restraints (kcal/mol)
    file: str
        name of .itp file for restraints (e.g. 'dihre.itp')
    """

    fc_rad_a = forceconst_A 
    fc_rad_b = forceconst_B 

    file = open(file, 'w')
    file.write('[ intermolecular_interactions ] \n[ bonds ] \n')
    file.write('; ai     aj    type     bA      fcA       bB      fcB\n')
    file.write(' %s   %s   6   %.2f   %.1f   %.2f   %.1f\n' % (
    atoms_restrain[0], atoms_restrain[1], dist, fc_rad_a, dist, fc_rad_b))

    file.close()

    return
