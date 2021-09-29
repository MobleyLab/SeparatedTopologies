###Calculate Boresch-style restraints in GROMACS
#Ligand atoms: largest ring systems, atoms closest to ligand COM
#protein atoms: middle of a helix, backbone/C-beta, >1nm away from ligand COM, check angles

import numpy as np
import mdtraj as md
import itertools
from simtk import unit
from scipy.spatial import distance
from openeye import oechem
import math
import sys

force_const = 83.68
R = 8.31445985*0.001  # Gas constant in kJ/mol/K
T = 298.15
RT = R*T

def select_ligand_atoms(lig, traj, ligand='LIG'):
    """Select three ligand atoms for Boresch-stylre restraints.
    Parameters
    ----------
    lig: str
        mol2 file of the ligand
    traj : str
        Mdtraj object with coordinates of the solvated protein-ligand system (e.g. from .gro file)
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

    #Load in ligand from mol2 into openeye
    ifs = oechem.oemolistream(lig)
    mol = oechem.OEGraphMol()
    oechem.OEReadMolecule(ifs, mol)

    # Find ring systems (here not specifically aromatic, could change that)
    nr_ring_systems, parts = oechem.OEDetermineRingSystems(mol)
    ring_systems = []
    atoms = []
    # If we find a ring system in our molecule
    if nr_ring_systems >= 1:
        #Loop through ringsystems, count number of atoms
        for ringidx in range(1, nr_ring_systems + 1):
            single_ring = []
            nr_atoms = parts.count(ringidx)
            ring_systems.append(nr_atoms)
            #store atom indices of ring systems
            for atom in mol.GetAtoms():
                if parts[atom.GetIdx()] == ringidx:
                    single_ring.append(atom.GetIdx())
            atoms.append(single_ring)
        #find largest ring/rings
        largest_ring = [i for i, x in enumerate(ring_systems) if x == max(ring_systems)]
        #Find indices of atoms in largest ring/rings
        ring_atoms = []
        for lr in largest_ring:
            ring_atoms.extend(atoms[lr])
        #select atoms of that ring system to use for restraints
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
    ###If there are multiple rings, could ligand atoms be in different rings??? Change that?
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
    sec_struc:  str
        secondary structure element: H for helix, E for strand
    residues2exclude: list
        list of residues (index) to exclude for Boresch atom selection (e.g. flexible loop)
        default = None
    Returns
    -------
    protein_list: list
        list of indices of possible protein atoms
                """
    topology = traj.topology

    if residues2exclude==None:
        heavy_protein = topology.select('protein and (backbone or name CB)')
        heavy_protein_traj = traj.atom_slice(heavy_protein, inplace=False)
        #Compute secondary structure of residues, output numpy array
        dssp = md.compute_dssp(heavy_protein_traj, simplified=True)
        structure = dssp[0].tolist()

        # Loop over all residues, look for start of a helix
        ex = []
        start_helix = False
        helix = []
        # How many residues of the protein to discard at the beginning and the end since ends can be floppy
        skip_start = 20
        skip_end = 10
        # number of residues Helix/beta sheet has to consist of to be considered stable
        stable_helix = 8

        #Check if more helices, more beta sheets in structure, choose predominant one
        res_in_helix = structure.count('H')
        res_in_sheet = structure.count('E')
        if res_in_helix >= res_in_sheet:
            sec_struc = ['H']
        else:
            sec_struc = ['H', 'E']
        for inx, b in enumerate(structure):

            if b in sec_struc and skip_start-1 < inx < (len(structure) - skip_end):
                # look for a start of a helix
                if structure[inx - 1] != b:

                    if structure[inx + 1:inx + stable_helix].count(b) == stable_helix-1:
                        start_helix = True
                        helix = []
                        helix.append('resid ' + str(inx))

                    else:
                        start_helix = False

                # Find end of helix
                elif structure[inx - 4:inx + 1].count(b) == 5 and structure[inx + 1] != b and start_helix == True:
                    helix.append('resid ' + str(inx))
                    # Leave out first 3 and last 3 residues of loop/sheet
                    ex.append(helix[3:-3])

                # If structure ends with helix account for that
                elif structure[inx - 4:inx + 1].count(b) == 5 and inx + 1 == (len(structure) - 6) and start_helix == True:
                    helix.append('resid ' + str(inx))
                    # Leave out first 3 and last 3 residues of loop/sheet
                    ex.append(helix[3:-3])
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

    #Select atoms that are at least 1nm away from the ligand
    heavy_protein_l1_pairs = np.array(list(itertools.product(heavy_protein, [l1])))
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

def check_angle(angle):
    # check if angle is <5kT from 0 or 180
    check1 = 0.5 * force_const * np.power((angle - 0.0) / 180.0 * np.pi, 2)
    check2 = 0.5 * force_const * np.power((angle - 180.0) / 180.0 * np.pi, 2)
    ang_check_1 = check1 / RT
    ang_check_2 = check2 / RT
    if ang_check_1 < 5.0 or ang_check_2  < 5.0:
        return False
    return True

def substructure_search(mol2_lig, smarts):
    """Pick ligand atoms for restraints based on substructure search.
        Parameters
        ----------
        mol2_lig : str
           mol2 file of the ligand.
        smarts : str
            Smarts pattern used for substructure search. One atom is flagged and picked for restraints.
        Returns
        -------
        matches: list
            Indices of ligand atoms that match the picked atom in the substructure.
        """
    ifs = oechem.oemolistream(mol2_lig)

    mol = oechem.OEGraphMol()

    oechem.OEReadMolecule(ifs, mol)

    query = oechem.OEQMol()
    oechem.OEParseSmarts(query, smarts)
    substructure_search = oechem.OESubSearch(query)

    matches = []

    for match in substructure_search.Match(mol, True):
        for matched_atom in match.GetAtoms():
            if matched_atom.pattern.GetMapIdx() != 0:
                # add index of matched atom/atoms to matches
                matches.append(matched_atom.target.GetIdx())

    return matches

def select_Boresch_atoms(traj, mol2_lig, ligand_atoms = None, protein_atoms = None, substructure = None, ligand='LIG'):
    """Select possible protein atoms for Boresch-style restraints.
    Parameters
    ----------
    traj : str
        Mdtraj object with coordinates of the system (e.g. from .gro file)
    mol2_lig: str
        mol2 file of the ligand
    ligand_atoms: list
        manual selection of three ligand atoms to choose for restraints, list of three indices of ligand atoms
    protein_atoms: list
        manual selection of three protein atoms to choose for restraints, list of three indices of protein atoms
        P1, P2, P3: P3 is used for distance restraint with L1
        Index of protein atoms is zero based
    substructure: list
        List of three strings of SMARTS pattern, each having one atom tagged. Atoms are tagged with :1
        Ex: ['[#6X3:1]:[#7X2]:[#6X3]', '[#6X3]:[#7X2:1]:[#6X3]', '[#6X3]:[#7X2]:[#6X3:1]']
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
            raise ValueError('Three ligand atoms need to be defined')
    # If a substructure (three SMARTS pattern with tagged atoms) is provided by the user
    elif substructure != None:
        topology = traj.topology
        ligand = topology.select('resname %s' % ligand).tolist()
        lig_length = len(ligand)
        ligand_atoms = []
        # For each ligand atom different pattern: List should have length of three
        if len(substructure) == 3:
            for pattern in substructure:
                match = substructure_search(mol2_lig, pattern)
                if len(match) == 0:
                    raise ValueError('No match was found')
                if len(match) > 1:
                    print('More than one match found. Picking first match.')
                    ligand_atoms.append(ligand[match[0]])
                else:
                    ligand_atoms.append(ligand[match[0]])
            l1, l2, l3 = ligand_atoms[0], ligand_atoms[1], ligand_atoms[2]
        else:
            raise ValueError('Three smarts pattern have to be provided.')
    else:
        #Get ligand atoms through automatic selection
        l1, l2, l3, lig_length = select_ligand_atoms(mol2_lig, traj, ligand=ligand)

    #Protein atoms: should be backbone/CB, part of a helix
    proteinlist = protein_list(traj, l1)

    #If protein atoms are specified by the user: Check that those are appropriate
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
                print(collinear, check_a1, check_a2)
                print(a1, a2)
                print('Manual selection not appropriate. Continuing with automatic selection for protein atoms')
                protein_atoms = None
            # Check if user specified protein atoms are among 'stable ones'
            # If not, still use them but let user know that these might not be stable
            elif protein_atoms[0] not in proteinlist and protein_atoms[1] not in proteinlist and protein_atoms[2] not in proteinlist:
                print('These protein atoms might not be a good selection. Check them (backbone?part of helix?).')

    # If no protein atoms are defined: Get protein atoms through automatic selection
    if protein_atoms == None:
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

        # for P2: choose first atom that passes checks and is at minimum distance of 5A from p1
        min_distance = 5 * unit.angstrom / unit.nanometer

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
        # to do: maybe just choose min distance from P1 and P2 instead of the one furthest away?
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

    values = [d[0][0], a1[0][0], a2[0][0], dih1[0][0], dih2[0][0], dih3[0][0]]
    restrained_atoms = [i + 1 for i in restrained_atoms]
    print('Boresch values:', values)
    print('Boresch atoms:', restrained_atoms)

    return values, restrained_atoms


def write_itp_restraints(restrained_atoms, values, fc_A, fc_B, file):
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

    file = open(file, 'w')
    file.write('[ intermolecular_interactions ] \n[ bonds ] \n')
    file.write('; ai     aj    type   bA      kA     bB      kB\n')
    file.write(' %s   %s   6   %.3f   %.1f   %.3f   %.1f\n\n' % (
    restrained_atoms[2], restrained_atoms[3], values[0], fc_A[0], values[0], fc_B[0]))
    file.write('[ angles ]\n')
    file.write('; ai     aj    ak     type    thA      fcA        thB      fcB\n')
    file.write(' %s   %s   %s   1   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[1], restrained_atoms[2], restrained_atoms[3], values[1], fc_A[1], values[1], fc_B[1]))
    file.write(' %s   %s   %s   1   %.2f   %.2f   %.2f   %.2f\n\n' % (
    restrained_atoms[2], restrained_atoms[3], restrained_atoms[4], values[2], fc_A[2], values[2], fc_B[2]))
    file.write('[ dihedrals ]\n')
    file.write('; ai     aj    ak    al    type     thA      fcA       thB      fcB\n')
    file.write(' %s   %s   %s   %s   2   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[0], restrained_atoms[1], restrained_atoms[2], restrained_atoms[3], values[3], fc_A[3], values[3], fc_B[3]))
    file.write(' %s   %s   %s   %s   2   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[1], restrained_atoms[2], restrained_atoms[3], restrained_atoms[4], values[4], fc_A[4], values[4], fc_B[3]))
    file.write(' %s   %s   %s   %s   2   %.2f   %.2f   %.2f   %.2f\n' % (
    restrained_atoms[2], restrained_atoms[3], restrained_atoms[4], restrained_atoms[5], values[5], fc_A[5], values[5], fc_B[3]))

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

def dist_correction_fc_angle(dist_x):
    # the length of the distance has an impact on the mobility of the ligand (arc length of sperical shell that
    # it can move on --> increase the force constant of the angle 1 quadratically with the distance
    # at a distance of 5A --> fc = 40 kcal/mol*rad2

    dist_0 = 0.5

    fc_0 = 167.36 #40 kcal/mol*rad2 at a distance of 5A

    fc_x = round((dist_x/dist_0)**2 * fc_0,2)

    return fc_x

def restrain_ligands(complex_A, complex_B, mol2_ligA, mol2_ligB, file_A0, file_B0, file_A1, file_B1,
                     ligand_atoms=None, protein_atoms=None, substructure = None, ligand='LIG'):
    """Select possible protein atoms for Boresch-style restraints, write .itp files for restraints.
     Parameters
     ----------
    complex_A : str
        coordinates of the system (e.g. from .gro file)
    complex_B : str
        coordinates of the system (e.g. from .gro file)
    mol2_ligA: str
        mol2 file of the ligand A
    mol2_ligB: str
        mol2 file of the ligand B
    file_A0: str
        name of .itp file to turn on restraints ligand A
    file_B0: str
        name of .itp file to turn off restraints ligand B
    file_A1: str
        name of .itp file for ligand A, restraints on in A and B state
    file_B1: str
        name of .itp file for ligand B, restraints on in A and B state
    ligand_atoms: list
        manual selection of three ligand atoms to choose for restraints, list of three indices of ligand atoms
    protein_atoms: list
        manual selection of three protein atoms to choose for restraints, list of three indices of protein atoms
    substructure: list
        List of three strings of SMARTS pattern, each having one atom tagged. Atoms are tagged with :1
        Ex: ['[#6X3:1]:[#7X2]:[#6X3]', '[#6X3]:[#7X2:1]:[#6X3]', '[#6X3]:[#7X2]:[#6X3:1]']
    ligand: str
        three letter code for the ligand
     Returns
     -------
     restrained_atoms : list
         List of indices of 3 protein and 3 ligand atoms selected for Boresch restraints for ligand A and ligand B
     """
    complex_A = md.load(complex_A)
    complex_B = md.load(complex_B)
    ###Restrained atoms
    ###do that for both ligands separately, then change the ligand atoms of ligandB since in the combined .gro they are different
    ###Indices for protein and ligand atoms are 0 based (index from file - 1)
    restrained_atoms_A, ligA_length = select_Boresch_atoms(complex_A, mol2_ligA, ligand_atoms=ligand_atoms, protein_atoms=protein_atoms, substructure=substructure, ligand=ligand)
    restrained_atoms_B, ligB_length = select_Boresch_atoms(complex_B, mol2_ligB, ligand_atoms=ligand_atoms, protein_atoms=restrained_atoms_A[:3], substructure=substructure, ligand=ligand)

    ###Compute distance, angles, dihedrals
    values_A, restrained_atoms_A = compute_dist_angle_dih(complex_A, restrained_atoms_A)
    values_B, restrained_atoms_B = compute_dist_angle_dih(complex_B, restrained_atoms_B)


    # For ligand B add length of ligand A since in combined .gro file
    restrained_atoms_B = edit_indices_ligandB(restrained_atoms_B, ligA_length)

    # Typically we restrain everything with 20 kcal/mol; here given in kJ/mol
    fc = 20 * 4.184

    # Get the force constant for the angle 1 which depends on the length of the distance restraint
    fc_thA_A = dist_correction_fc_angle(values_A[0])
    fc_thA_B = dist_correction_fc_angle(values_B[0])

    dG_A_off_analytical = analytical_Boresch_correction(values_A[0], values_A[1], values_A[2], fc, fc_thA_A,
                                                        fc, fc, fc, fc, T=298.15)

    dG_B_off_analytical = analytical_Boresch_correction(values_B[0], values_B[1], values_B[2], fc, fc_thA_B,
                                                        fc, fc, fc, fc, T=298.15)
    dG_B_on_analytical = -dG_B_off_analytical
    ###write .itp for restraints section
    #Restraining: A state fc_A = 0
    #Restrain distance, angle and dihedral with a force constant of 20 kcal/mol*A2 / kcal/mol*rad2

    #Force constant for distance, 2 angles, 3 dihedrals
    force_const_lig_A = [fc*100, fc_thA_A, fc, fc, fc, fc]
    force_const_lig_B = [fc*100, fc_thA_B, fc, fc, fc, fc]
    force_const_0 = [0,0,0,0,0,0]
    write_itp_restraints(restrained_atoms_A, values_A, force_const_0, force_const_lig_A, file_A0)
    write_itp_restraints(restrained_atoms_B, values_B, force_const_lig_B, force_const_0, file_B0)


    #FEC: fc_A = fc_B
    write_itp_restraints(restrained_atoms_A, values_A, force_const_lig_A, force_const_lig_A, file_A1)
    write_itp_restraints(restrained_atoms_B, values_B, force_const_lig_B, force_const_lig_B, file_B1)
    return restrained_atoms_A, restrained_atoms_B, round(dG_A_off_analytical,3), round(dG_B_on_analytical,3)


def analytical_Boresch_correction(r0, thA, thB, fc_r, fc_thA, fc_thB, fc_phiA, fc_phiB, fc_phiC, T=298.15):
    #Analytical correction orientational restraints, equation 14 from Boresch 2003 paper doi 10.1021/jp0217839
    K = 8.314472 * 0.001  # Gas constant in kJ/mol/K
    V = 1.66  # standard volume in nm^3
    # distance: Force constants in kJ/mol*nm2
    fc_r = fc_r * 100

    thA = math.radians(thA)  # get angle in radians
    thB = math.radians(thB)  # get angle in radians

    # dG for turning Boresch restraints off (switch sign if you want to turn them on)
    dG = - K * T * math.log(((8.0 * math.pi ** 2.0 * V) / (r0 ** 2.0 * math.sin(thA) * math.sin(thB))
            *
            (((fc_r * fc_thA * fc_thB * fc_phiA * fc_phiB * fc_phiC) ** 0.5) / ((2.0 * math.pi * K * T) ** (3.0)))))

    #dG in kcal mol
    dG = dG / 4.184
    return dG