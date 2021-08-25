###Calculate Boresch-style restraints in GROMACS
#ligand atoms: heavy atoms, L1 closest to COM of ligand, L2,L3 close to L1
#protein atoms: backbone and CB, R3 within 1-5A of L1, R1 and R2 within 1-8A
#check that atoms not collinear


import numpy as np
import mdtraj as md
import itertools
from simtk import unit
from scipy.spatial import distance
import random

random.seed(8)

def select_ligand_atoms(traj, ligand='LIG'):
    """Select three ligand atoms for Boresch-stylre restraints.
    Parameters
    ----------
    traj : mdtraj trajectory
        Mdtraj object with coordinates of the system (e.g. from .gro file)
    ligand : str
        Three letter code for ligand name
    Returns
    -------
    L1, L2, L3 : int
        Index of ligand atoms selected for Boresch restraints
    lig_length: int
        Number of atoms of the ligand
    """
    topology = traj.topology
    #Only consider heavy atoms for restraints
    heavy_ligand = topology.select('resname %s and not element H'%ligand).tolist()
    ligand = topology.select('resname %s'%ligand).tolist()
    lig_length = len(ligand)
    ligand_traj = traj.atom_slice(ligand, inplace=False)
    #Compute the center of mass of the ligand
    com_ligand = md.compute_center_of_mass(ligand_traj)
    ligand_coordinates = traj.xyz[:, heavy_ligand, :]
    com_distance = []
    #Calculate the distance of ligand atoms to the COM
    for c in ligand_coordinates[0]:
        d = distance.euclidean(com_ligand[0], c)
        com_distance.append(d)
    # store index with sorted value
    sorted_distance = sorted((x, y) for y, x in enumerate(com_distance))

    #Find ligand heavy atom closest to COM
    index_L1 = sorted_distance[0][1]
    coordinates_L1 = ligand_coordinates[0][index_L1]
    L1 = heavy_ligand[index_L1]

    ###Find 2 other heavy atoms near first carbon in ligand
    L1_distance = []
    for c in ligand_coordinates[0]:
        d = distance.euclidean(coordinates_L1, c)
        L1_distance.append(d)
    sorted_distance = sorted((x, y) for y, x in enumerate(L1_distance))

    index_L2 = sorted_distance[1][1]
    index_L3 = sorted_distance[2][1]
    L2 = heavy_ligand[index_L2]
    L3 = heavy_ligand[index_L3]

    return L1, L2, L3, lig_length

def r3_r12_list(traj, l1, residues2exclude=None):
    """Select possible protein atoms for Boresch-stylre restraints.
    Parameters
    ----------
    traj : mdtraj trajectory
        Mdtraj object with coordinates of the system (e.g. from .gro file)
    l1 : int
        Index of ligand atom L1 chosen for restraints
    residues2exclude: list
        list of residues (index) to exclude for Boresch atom selection
        default = None
    Returns
    -------
    r3 : int
        Index of protein atom R3 selected for Boresch restraints
    r12_list: list
        list of indices of possible R1 and R2 protein atoms
                """
    topology = traj.topology
    #select backbone and CB atoms protein
    if residues2exclude==None:
        heavy_protein = topology.select('protein and (backbone or name CB)').tolist()

    #if a list of residue indices is provided: exclude those residues
    else:
        ex = []
        for index, v in enumerate(residues2exclude):
            exclude = 'resid ' + str(v)
            if index < len(residues2exclude)-1:
                exclude = exclude + ' or '
            ex.append(exclude)
        ex = ''.join(ex)
        heavy_protein = topology.select('protein and (backbone or name CB) and not (%s)'%ex).tolist()
    ###Find protein atom R3 within distance of L1
    # next couple of lines modified from yank.restraints.py
    heavy_protein_l1_pairs = np.array(list(itertools.product(heavy_protein, [l1])))

    # Filter r3-l1 pairs that are too close/far away for the distance constraint.
    max_distance_r12 = 8 * unit.angstrom / unit.nanometer
    max_distance = 5 * unit.angstrom / unit.nanometer
    min_distance = 1 * unit.angstrom / unit.nanometer

    distances = md.geometry.compute_distances(traj, heavy_protein_l1_pairs)[0]

    indices_of_in_range_pairs = np.where(np.logical_and(distances > min_distance, distances <= max_distance))[0]
    r3_l1_pairs = heavy_protein_l1_pairs[indices_of_in_range_pairs].tolist()
    #Choose a random r3 protein atom that is within 1 - 5 A from L1
    r3 = random.sample(r3_l1_pairs, 1)[0][0]

    #R1 and R2 protein atom that is within 1 - 8 A from L1
    indices_of_in_range_pairs_r12 = np.where(np.logical_and(distances > min_distance, distances <= max_distance_r12))[0]

    r12_list = [heavy_protein[i] for i in indices_of_in_range_pairs_r12]

    return r3, r12_list

#copied from yank version 0.25.2
def _is_collinear(positions, atoms, threshold=0.9):
    """Report whether any sequential vectors in a sequence of atoms are collinear.
    Parameters
    ----------
    positions : n_atoms x 3 simtk.unit.Quantity
        Reference positions to use for imposing restraints (units of length).
    atoms : iterable of int
        The indices of the atoms to test.
    threshold : float, optional, default=0.9
        Atoms are not collinear if their sequential vector separation dot
        products are less than ``threshold``.
    Returns
    -------
    result : bool
        Returns True if any sequential pair of vectors is collinear; False otherwise.
    """
    result = False
    for i in range(len(atoms) - 2):
        v1 = positions[atoms[i + 1], :] - positions[atoms[i], :]
        v2 = positions[atoms[i + 2], :] - positions[atoms[i + 1], :]
        normalized_inner_product = np.dot(v1, v2) / np.sqrt(np.dot(v1, v1) * np.dot(v2, v2))
        result = result or (normalized_inner_product > threshold)

    return result

def atoms_2_restrain(traj, ligand='LIG'):
    """Select possible protein atoms for Boresch-stylre restraints.
    Parameters
    ----------
    traj : mdtraj trajectory
        Mdtraj object with coordinates of the system (e.g. from .gro file)
    Returns
    -------
    restrained_atoms : list
        List of indices of 3 protein and 3 ligand atoms selected for Boresch restraints
    lig_length: int
        Number of ligand atoms
    """
    #Get ligand atoms
    l1, l2, l3, lig_length = select_ligand_atoms(traj, ligand=ligand)
    #Get protein atoms
    r3, r12_list = r3_r12_list(traj, l1)
    complex_coordinates = traj.xyz[:, :, :]
    ###Iterate over different R1 and R2 atoms and check collinearity
    restrained_atoms_list = []
    for r1 in r12_list:

        for r2 in r12_list:

            restrained_atoms = [r1, r2, r3, l1, l2, l3]
            if r1 != r2 and r1 != r3 and r2 != r3:
                collinear = _is_collinear(complex_coordinates[0], restrained_atoms)
                #only choose atoms that are not collinear
                if collinear == False:
                    restrained_atoms_list.append(restrained_atoms)

                else:
                    continue

            else:
                continue

    restrained_atoms = random.sample(restrained_atoms_list, 1)[0]
    # Add one since python starts at 0, .gro file with 1
    return restrained_atoms, lig_length

def compute_dist_angle_dih(complex, restrained_atoms):
    """Compute distance, angles, dihedrals.
    Parameters
    ----------
    complex : mdtraj trajectory
        Mdtraj object with coordinates of the system (e.g. from .gro file)
    restrained_atoms: list
        list of protein and ligand atoms selected for restraints
    Returns
    -------
    values : list
        List of values for distance, two angles and three dihedral
    """
    d = md.compute_distances(complex, [[restrained_atoms[2], restrained_atoms[3]]])
    a1 = np.rad2deg(md.geometry.compute_angles(complex, np.array([restrained_atoms[1:4]])))
    a2 = np.rad2deg(md.geometry.compute_angles(complex, np.array([restrained_atoms[2:5]])))
    dih1 = np.rad2deg(md.compute_dihedrals(complex, [np.array(restrained_atoms[0:4])]))
    dih2 = np.rad2deg(md.compute_dihedrals(complex, [np.array(restrained_atoms[1:5])]))
    dih3 = np.rad2deg(md.compute_dihedrals(complex, [np.array(restrained_atoms[2:6])]))

    values = [d, a1, a2, dih1, dih2, dih3]
    # print(values)
    restrained_atoms = [i + 1 for i in restrained_atoms]
    # print(restrained_atoms)
    return values, restrained_atoms


def write_itp_restraints(restrained_atoms, values, forceconst_A, forceconst_B, file):
    """Compute distance, angles, dihedrals.
    Parameters
    ----------
    restrained_atoms: list
        list of protein and ligand atoms selected for restraints
    values: list
        List of values for distance, two angles and three dihedral
    forceconst: int
        forceconstant for restraints
    file: str
        name of .itp file for restraints (e.g. 'boresch_restraints_A.itp' for ligand A)
    fc_A: bool
        True: A-state is restrained
        False: A-state is unrestrained (fc=0), B-state is restrained
    """

    fc_dist_a = forceconst_A * 100 * 4.184
    fc_rad_a = forceconst_A * 4.184
    fc_dist_b = forceconst_B * 100 * 4.184
    fc_rad_b = forceconst_B * 4.184

    file = open(file, 'w')
    file.write('[ intermolecular_interactions ] \n[ bonds ] \n')
    file.write('; ai     aj    type   bA      kA     bB      kB\n')
    file.write(' %s   %s   6   %.3f   %.1f   %.3f   %.1f\n\n' % (
    restrained_atoms[2], restrained_atoms[3]+2, values[0], fc_dist_a, values[0], fc_dist_b))
    file.write('[ angles ]\n')
    file.write('; ai     aj    ak     type    thA      fcA        thB      fcB\n')
    file.write(' %s   %s   %s   1   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[1], restrained_atoms[2], restrained_atoms[3]+2, values[1], fc_rad_a, values[1], fc_rad_b))
    file.write(' %s   %s   %s   1   %.2f   %.2f   %.2f   %.2f\n\n' % (
    restrained_atoms[2], restrained_atoms[3]+2, restrained_atoms[4]+2, values[2], fc_rad_a, values[2], fc_rad_b))
    file.write('[ dihedrals ]\n')
    file.write('; ai     aj    ak    al    type     thA      fcA       thB      fcB\n')
    file.write(' %s   %s   %s   %s   2   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[0], restrained_atoms[1], restrained_atoms[2], restrained_atoms[3]+2, values[3], fc_rad_a, values[3], fc_rad_b))
    file.write(' %s   %s   %s   %s   2   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[1], restrained_atoms[2], restrained_atoms[3]+2, restrained_atoms[4]+2, values[4], fc_rad_a, values[4], fc_rad_b))
    file.write(' %s   %s   %s   %s   2   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[2], restrained_atoms[3]+2, restrained_atoms[4]+2, restrained_atoms[5]+2, values[5], fc_rad_a, values[5], fc_rad_b))

    file.close()

    return

def edit_indices_ligandB(restrained_atoms, ligA_length):
    # Indices of ligand B atoms change when transferring it to complex A --> account for that
    restrained_atoms = restrained_atoms[0:3] + [i + ligA_length for i in restrained_atoms[3:6]]
    return restrained_atoms

def include_itp_in_top(top, idpfile):
    with open(top) as file:
        newline = ''
        for line in file:
            for part in line.split():
                if 'boresch_restraints' not in part:
                    newline = '\n#include "%s"'%idpfile
        file.close()
    file = open(top, 'a')
    file.write(newline)
    file.close()

def restrain_ligands(complex_A, complex_B, file_A0, file_B0, file_A1, file_B1, ligand='LIG'):
    complex_A = md.load(complex_A)
    complex_B = md.load(complex_B)
    ###Restrained atoms
    ###do that for both ligands separately, then change the ligand atoms of ligandB since in the combined .gro they are different
    restrained_atoms_A, ligA_length = atoms_2_restrain(complex_A, ligand=ligand)
    restrained_atoms_B, ligB_length = atoms_2_restrain(complex_B, ligand=ligand)

    ###Compute distance, angles, dihedrals
    values_A, restrained_atoms_A = compute_dist_angle_dih(complex_A, restrained_atoms_A)
    values_B, restrained_atoms_B = compute_dist_angle_dih(complex_B, restrained_atoms_B)

    # For ligand B add length of ligand A since in combined .gro file
    restrained_atoms_B = edit_indices_ligandB(restrained_atoms_B, ligA_length)
    ###write .itp for restraints section
    #Restraining: A state fc_A = 0

    write_itp_restraints(restrained_atoms_A, values_A, 0, 20, file_A0)
    write_itp_restraints(restrained_atoms_B, values_B, 20, 0, file_B0)


    #FEC: fc_A = fc_B
    write_itp_restraints(restrained_atoms_A, values_A, 20, 20, file_A1)
    write_itp_restraints(restrained_atoms_B, values_B, 20, 20, file_B1)
    return restrained_atoms_A, restrained_atoms_B

def restrain_ligands_boresch_endstate(complex_A, complex_B, file_A0, file_B0, ligand='LIG'):
    complex_A = md.load(complex_A)
    complex_B = md.load(complex_B)
    ###Restrained atoms
    ###do that for both ligands separately, then change the ligand atoms of ligandB since in the combined .gro they are different
    restrained_atoms_A, ligA_length = atoms_2_restrain(complex_A, ligand=ligand)
    restrained_atoms_B, ligB_length = atoms_2_restrain(complex_B, ligand=ligand)

    ###Compute distance, angles, dihedrals
    values_A, restrained_atoms_A = compute_dist_angle_dih(complex_A, restrained_atoms_A)
    values_B, restrained_atoms_B = compute_dist_angle_dih(complex_B, restrained_atoms_B)

    ###write .itp for restraints section
    #Restraining: A state fc_A = 0
    # For ligand B add length of ligand A since in combined .gro file
    restrained_atoms_B = edit_indices_ligandB(restrained_atoms_B, ligA_length)
    write_itp_restraints(restrained_atoms_A, values_A, 0, 20, file_A0)
    write_itp_restraints(restrained_atoms_B, values_B, 20, 0, file_B0)

    return restrained_atoms_A, restrained_atoms_B