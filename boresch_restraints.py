###Calculate Boresch-style restraints in GROMACS
#ligand atoms: heavy atoms, L1 closest to COM of ligand, L2,L3 close to L1
#protein atoms: backbone and CB, R3 within 1-5A of L1, R1 and R2 within 1-8A
#check that atoms not collinear


import numpy as np
import mdtraj as md
import itertools
from simtk import unit
from scipy.spatial import distance
from openeye import oechem
import tempfile
import os

def select_ligand_atoms(traj, ligand='LIG'):
    """Select three ligand atoms for Boresch-stylre restraints.
    Parameters
    ----------
    traj : str
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
    # Only consider heavy atoms for restraints
    heavy_ligand = topology.select('resname %s and not element H' % ligand).tolist()
    ligand = topology.select('resname %s' % ligand).tolist()
    lig_length = len(ligand)
    ligand_traj = traj.atom_slice(ligand, inplace=False)

    # Get openeye molecule from mdtraj via pdb file (code adapted from perses)
    # create a temporary file with a PDB suffix and save with MDTraj
    pdb_file = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    ligand_traj.save(pdb_file.name)
    pdb_file.close()

    # Now use the openeye oemolistream to read in this file as an OEMol:
    ifs = oechem.oemolistream()
    ifs.open(pdb_file.name)
    ifs.SetFormat(oechem.OEFormat_PDB)
    lig = oechem.OEGraphMol()
    oechem.OEReadMolecule(ifs, lig)

    # close the stream and delete the temporary pdb file
    ifs.close()
    os.unlink(pdb_file.name)

    nr_ring_systems, parts = oechem.OEDetermineRingSystems(lig)
    ring_systems = []
    atoms = []
    if nr_ring_systems >= 1:
        for ringidx in range(1, nr_ring_systems + 1):
            single_ring = []
            nr_atoms = parts.count(ringidx)
            ring_systems.append(nr_atoms)
            #store atom indices of ring systems
            for atom in lig.GetAtoms():
                if parts[atom.GetIdx()] == ringidx:

                    single_ring.append(atom.GetIdx())
            atoms.append(single_ring)
        #find largest ring/rings
        largest_ring = [i for i, x in enumerate(ring_systems) if x == max(ring_systems)]
        #Find indices of atoms in largest ring/rings
        ring_atoms = []
        for lr in largest_ring:
            ring_atoms.extend(atoms[lr])

        ligand_atoms = [heavy_ligand[i] for i in ring_atoms]
        ligand_coordinates = traj.xyz[:, ligand_atoms, :]

    #If no ring systems in molecule just choose all heavy atoms
    elif nr_ring_systems == 0:
        ligand_coordinates = traj.xyz[:, heavy_ligand, :]
        ligand_atoms = heavy_ligand

    # Compute the center of mass of the ligand
    com_ligand = md.compute_center_of_mass(ligand_traj)

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
    L1 = ligand_atoms[index_L1]

    ###Find 2 other heavy atoms near first selected ligand atom
    L1_distance = []
    for c in ligand_coordinates[0]:
        d = distance.euclidean(coordinates_L1, c)
        L1_distance.append(d)
    sorted_distance = sorted((x, y) for y, x in enumerate(L1_distance))

    index_L2 = sorted_distance[1][1]
    index_L3 = sorted_distance[2][1]
    L2 = ligand_atoms[index_L2]
    L3 = ligand_atoms[index_L3]

    return L1, L2, L3, lig_length

def protein_list(traj, l1, residues2exclude=None):
    """Select possible protein atoms for Boresch-style restraints.
    Parameters
    ----------
    traj : str
        Mdtraj object with coordinates of the system (e.g. from .gro file)
    l1 : int
        Index of ligand atom L1 chosen for restraints
    residues2exclude: list
        list of residues (index) to exclude for Boresch atom selection (e.g. flexible loop)
        default = None
    Returns
    -------
    protein_list: list
        list of indices of possible protein atoms
                """
    topology = traj.topology
    #select backbone and CB atoms protein
    if residues2exclude==None:
        heavy_protein = topology.select('protein and (backbone or name CB)')
        heavy_protein_traj = traj.atom_slice(heavy_protein, inplace=False)
        #Compute secondary structure of residues, output numpy array
        dssp = md.compute_dssp(heavy_protein_traj)
        structure = dssp[0].tolist()

        # to DO: maybe also allow beta strands: in general of if no appropriate atoms are found using just helices
        # Loop over all residues, look for start of a helix
        ex = []
        start_helix = False
        helix = []
        for inx, b in enumerate(structure):

            # discard first 6 and last 6 residues of the protein since they can be floppy

            if b == 'H' and 6 < inx < (len(structure) - 6):
                # look for a start of a helix
                if structure[inx - 1] != 'H':
                    # helix must have at least 6 residues to be considered stable
                    if structure[inx + 1:inx + 6].count('H') == 5:
                        start_helix = True
                        helix = []
                        helix.append('resid ' + str(inx))


                    else:
                        start_helix = False

                # Find end of helix
                elif structure[inx - 4:inx + 1].count('H') == 5 and structure[inx + 1] != 'H':
                    end_helix = True
                    helix.append('resid ' + str(inx))
                    #Find middle of the helix
                    middle_helix = (int(len(helix)/2))
                    #Only choose middle residue + 2 surrounding
                    ex.append(helix[middle_helix - 3:middle_helix+2])

                # If structure ends with helix account for that
                elif structure[inx - 4:inx + 1].count('H') == 5 and inx + 1 == (len(structure) - 6):
                    helix.append('resid ' + str(inx))
                    # Find middle of the helix
                    middle_helix = (int(len(helix) / 2))
                    # Only choose middle residue + 2 surrounding
                    ex.append(helix[middle_helix - 3:middle_helix + 2])
                else:
                    if start_helix == True:
                        helix.append('resid ' + str(inx))

                    else:
                        continue
        ex = list(itertools.chain.from_iterable(ex))
        ex = ' '.join(ex)

        heavy_protein = topology.select('protein and (backbone or name CB) and (%s)' % ex).tolist()

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
    min_distance = 10 * unit.angstrom / unit.nanometer
    distances = md.geometry.compute_distances(traj, heavy_protein_l1_pairs)[0]

    indices_of_in_range_pairs = np.where(distances > min_distance)[0]

    protein_list = [heavy_protein[i] for i in indices_of_in_range_pairs]

    return protein_list

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

force_const = 83.68
R = 8.31445985*0.001  # Gas constant in kJ/mol/K
T = 298.15
RT = R*T

def check_angle(angle):
    # check if angle is <5kT from 0 or 180
    check1 = 0.5 * force_const * np.power((angle - 0.0) / 180.0 * np.pi, 2)
    check2 = 0.5 * force_const * np.power((angle - 180.0) / 180.0 * np.pi, 2)
    ang_check_1 = check1 / RT
    ang_check_2 = check2 / RT
    if ang_check_1 < 5.0 or ang_check_2  < 5.0:
        return False
    return True

def select_Boresch_atoms(traj, ligand_atoms = None, protein_atoms = None, ligand='LIG'):
    """Select possible protein atoms for Boresch-style restraints.
    Parameters
    ----------
    traj : str
        Mdtraj object with coordinates of the system (e.g. from .gro file)
    ligand_atoms: list
        manual selection of three ligand atoms to choose for restraints, list of three indices of ligand atoms
    protein_atoms: list
        manual selection of three protein atoms to choose for restraints, list of three indices of protein atoms
    ligand: str
        three letter code for the ligand
    Returns
    -------
    restrained_atoms : list
        List of indices of 3 protein and 3 ligand atoms selected for Boresch restraints
    lig_length: int
        Number of ligand atoms
    """

    complex_coordinates = traj.xyz[:, :, :]

    #If ligand atoms are user defined
    if ligand_atoms != None:
        if len(ligand_atoms) == 3:
            topology = traj.topology
            ligand = topology.select('resname %s' % ligand).tolist()
            lig_length = len(ligand)
            l1, l2, l3 = ligand_atoms[0], ligand_atoms[1], ligand_atoms[2]
        else:
            # Three ligand atoms are needed to define restraints
            print('Three ligand atoms need to be defined')
    else:
        #Get ligand atoms through automatic selection
        l1, l2, l3, lig_length = select_ligand_atoms(traj, ligand=ligand)

    #Protein atoms: should be backbone/CB, part of a helix
    proteinlist = protein_list(traj, l1)
    if protein_atoms != None:
        if len(protein_atoms) == 3:
            restrained_atoms = [protein_atoms[0], protein_atoms[1], protein_atoms[2], l1, l2, l3]
            #Make sure that user defined protein atoms pass collinearity criteria and angle checks
            collinear = _is_collinear(complex_coordinates[0], restrained_atoms)
            a1 = np.rad2deg(md.geometry.compute_angles(traj, np.array([[protein_atoms[1], protein_atoms[2], l1]])))
            check_a1 = check_angle(a1)
            a2 = np.rad2deg(md.geometry.compute_angles(traj, np.array([[protein_atoms[2],l1,l2]])))
            check_a2 = check_angle(a2)
            if collinear == True or check_a1 == False or check_a2 == False:
                print('Manual selection not appropriate. Continuing with automatic selection for protein atoms')
                protein_atoms = None
            # Check if user specified protein atoms are among 'stable ones'
            # If not, still use them but let user know that these might not be stable
            elif protein_atoms[0] not in proteinlist and protein_atoms[1] not in proteinlist and protein_atoms[2] not in proteinlist:
                print('These protein atoms might not be a good selection. Check them (backbone?part of helix?).')


    if protein_atoms == None:
    #Get protein atoms through automatic selection

        # pick P1 from this list that passes check_angle
        # choose first atom that we find
        # p1s = []
        for p in proteinlist:
            angle_1 = [p, l1, l2]
            a1 = np.rad2deg(md.geometry.compute_angles(traj, np.array([angle_1])))
            check_a1 = check_angle(a1)

            #check for collinearlity as well (make sure dihedral not 180)
            dih_1 = [p, l1, l2, l3]
            collinear = _is_collinear(complex_coordinates[0], dih_1)
            if check_a1 == True and collinear == False:
                p1 = p
                break


        #         p1s.append(p)
        # #pick P1 that is furthest from L1
        # dist_p1_l1 = []
        # for p1 in p1s:
        #     d = md.compute_distances(traj, [[p1, l1]])
        #     dist_p1_l1.append(d)
        # sorted_distance = sorted((x, y) for y, x in enumerate(dist_p1_l1))
        # index_P1 = sorted_distance[-1][1]
        # p1 = p1s[index_P1]

        # for P2: loop through list and save everything that passes check angle P2-P1-L1
        # choose first atom that's at minimum distance of 5A from p1
        min_distance = 5 * unit.angstrom / unit.nanometer


        # p2s = []
        for p in proteinlist:
            angle_2 = [p, p1, l1]
            a2 = np.rad2deg(md.geometry.compute_angles(traj, np.array([angle_2])))
            check_a2 = check_angle(a2)
            #check for collinearlity as well (make sure dihedral not 180)
            dih_2 = [p, p1, l1, l2, l3]
            collinear = _is_collinear(complex_coordinates[0], dih_2)
            if check_a2 == True and collinear == False:
                distance = md.geometry.compute_distances(traj, np.array([[p,p1]]))[0]

                if distance > min_distance:
                    p2 = p
                    break


                # p2s.append(p)
        # choose atom furthest from P1
        # dist_p1_p2 = []
        # for p2 in p2s:
        #     d = md.compute_distances(traj, [[p1, p2]])
        #     dist_p1_p2.append(d)
        # sorted_distance = sorted((x, y) for y, x in enumerate(dist_p1_p2))
        #
        # index_P2 = sorted_distance[-1][1]
        # p2 = p2s[index_P2]

        # for P3:loop through list and save everything that passes collinearity test
        p3s = []
        for p in proteinlist:
            restrained_atoms = [p, p2, p1, l1, l2, l3]

            if p != p2 and p != p1:
                collinear = _is_collinear(complex_coordinates[0], restrained_atoms)
                if collinear == False:
                    p3s.append(p)
                else:
                    continue
            else:
                continue
        # choose atom furthest from P1 and P2
        # to do: maybe just choose min distance from P1 and P2?
        dist_p1_p3 = []
        dist_p2_p3 = []
        for p in p3s:
            d = md.compute_distances(traj, [[p1, p]])
            dist_p1_p3.append(d)
            d = md.compute_distances(traj, [[p2, p]])
            dist_p2_p3.append(d)
        #the maximum of the product of the distances should be where the atom farthest from P1 and P2 is
        products = [a * b for a, b in zip(dist_p1_p3, dist_p2_p3)]
        sorted_distance = sorted((x, y) for y, x in enumerate(products))
        index_P3 = sorted_distance[-1][1]
        p3 = p3s[index_P3]


        restrained_atoms = [p3,p2,p1,l1,l2,l3]
        # Add one since python starts at 0, .gro file with 1
    return restrained_atoms, lig_length


# def atoms_2_restrain(traj, ligand='LIG'):
#     """Select possible protein atoms for Boresch-style restraints.
#     Parameters
#     ----------
#     traj : str
#         Mdtraj object with coordinates of the system (e.g. from .gro file)
#     Returns
#     -------
#     restrained_atoms : list
#         List of indices of 3 protein and 3 ligand atoms selected for Boresch restraints
#     lig_length: int
#         Number of ligand atoms
#     """
#     #Get ligand atoms
#     l1, l2, l3, lig_length = select_ligand_atoms(traj, ligand=ligand)
#     #Get protein atoms
#     r3, r12_list = r3_r12_list(traj, l1)
#     complex_coordinates = traj.xyz[:, :, :]
#     ###Iterate over different R1 and R2 atoms and check collinearity
#     restrained_atoms_list = []
#     for r1 in r12_list:
#
#         for r2 in r12_list:
#
#             restrained_atoms = [r1, r2, r3, l1, l2, l3]
#             if r1 != r2 and r1 != r3 and r2 != r3:
#                 collinear = _is_collinear(complex_coordinates[0], restrained_atoms)
#                 a1 = np.rad2deg(md.geometry.compute_angles(traj, np.array([restrained_atoms[1:4]])))
#                 a2 = np.rad2deg(md.geometry.compute_angles(traj, np.array([restrained_atoms[2:5]])))
#
#                 check_a1 = check_angle(a1)
#                 check_a2 = check_angle(a2)
#
#                 #only choose atoms that are not collinear and >=5kT from 0 or 180
#                 if collinear == False and check_a1 == True and check_a2 == True:
#                     restrained_atoms_list.append(restrained_atoms)
#
#                 else:
#                     continue
#
#             else:
#                 continue
#
#     restrained_atoms = random.sample(restrained_atoms_list, 1)[0]
#     # Add one since python starts at 0, .gro file with 1
#     return restrained_atoms, lig_length

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
    print(values)
    restrained_atoms = [i + 1 for i in restrained_atoms]
    print(restrained_atoms)
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
    restrained_atoms[2], restrained_atoms[3], values[0], fc_dist_a, values[0], fc_dist_b))
    file.write('[ angles ]\n')
    file.write('; ai     aj    ak     type    thA      fcA        thB      fcB\n')
    file.write(' %s   %s   %s   1   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[1], restrained_atoms[2], restrained_atoms[3], values[1], fc_rad_a, values[1], fc_rad_b))
    file.write(' %s   %s   %s   1   %.2f   %.2f   %.2f   %.2f\n\n' % (
    restrained_atoms[2], restrained_atoms[3], restrained_atoms[4], values[2], fc_rad_a, values[2], fc_rad_b))
    file.write('[ dihedrals ]\n')
    file.write('; ai     aj    ak    al    type     thA      fcA       thB      fcB\n')
    file.write(' %s   %s   %s   %s   2   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[0], restrained_atoms[1], restrained_atoms[2], restrained_atoms[3], values[3], fc_rad_a, values[3], fc_rad_b))
    file.write(' %s   %s   %s   %s   2   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[1], restrained_atoms[2], restrained_atoms[3], restrained_atoms[4], values[4], fc_rad_a, values[4], fc_rad_b))
    file.write(' %s   %s   %s   %s   2   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[2], restrained_atoms[3], restrained_atoms[4], restrained_atoms[5], values[5], fc_rad_a, values[5], fc_rad_b))

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

def restrain_ligands(complex_A, complex_B, file_A0, file_B0, file_A1, file_B1, ligand_atoms=None, protein_atoms=None, ligand='LIG'):
    complex_A = md.load(complex_A)
    complex_B = md.load(complex_B)
    ###Restrained atoms
    ###do that for both ligands separately, then change the ligand atoms of ligandB since in the combined .gro they are different
    ###Indices for protein and ligand atoms are 0 based (index from file - 1)
    restrained_atoms_A, ligA_length = select_Boresch_atoms(complex_A, ligand_atoms=ligand_atoms, protein_atoms=protein_atoms, ligand=ligand)
    restrained_atoms_B, ligB_length = select_Boresch_atoms(complex_B, ligand_atoms=ligand_atoms, protein_atoms=restrained_atoms_A[:3], ligand=ligand)

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