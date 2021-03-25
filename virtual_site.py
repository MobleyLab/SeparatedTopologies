import numpy as np
import mdtraj as md
import itertools
import shutil
from simtk import unit
from scipy.spatial import distance
import random

random.seed(8)


def _is_collinear(atoms, threshold=0.98):
    """Report whether any sequential vectors in a sequence of atoms are collinear.
    Parameters
    ----------
    atoms : list
        The 3D coordinates of the atoms to test.
    threshold : float, optional, default=0.98
        Atoms are not collinear if their sequential vector separation dot
        products are less than ``threshold``.
    Returns
    -------
    result : bool
        Returns True if any sequential pair of vectors is collinear; False otherwise.
    """
    result = False
    for i in range(len(atoms) - 2):
        v1 = atoms[i + 1] - atoms[i]
        v2 = atoms[i + 2] - atoms[i + 1]
        # v1 = positions[atoms[i + 1], :] - positions[atoms[i], :]
        # v2 = positions[atoms[i + 2], :] - positions[atoms[i + 1], :]
        normalized_inner_product = np.dot(v1, v2) / np.sqrt(np.dot(v1, v1) * np.dot(v2, v2))
        result = result or (normalized_inner_product > threshold)

    return result


def P1_virtual_site(protein_coordinates, heavy_protein, com_ligand):
    """Find first protein atom that is used to define the VS.
    Parameters
    ----------
    protein_coordinates : list
    heavy_protein : list
        Indices of heavy protein atoms.
    com_ligand: list
        Coordinates center of mass ligand
    Returns
    -------
    coordinates_P1
    P1+1: index protein atom
    distance_P1_COM: distance between protein atom and COM ligand
    """

    com_distance = []
    # Calculate the distance of protein atoms to the COM
    for c in protein_coordinates:
        d = distance.euclidean(com_ligand, c)
        com_distance.append(d)
        # store index with sorted value
    sorted_distance = sorted((x, y) for y, x in enumerate(com_distance))
    # For first reference protein atom for VS choose protein atom that is closest to COM ligand
    index_P1 = sorted_distance[0][1]
    distance_P1_COM = sorted_distance[0][0]
    coordinates_P1 = protein_coordinates[index_P1]
    P1 = heavy_protein[index_P1]

    return coordinates_P1, P1 + 1, distance_P1_COM


def P2_virtual_site(P2_list, complex_coordinates, com_ligand):
    """ Find second protein atom that is used to define the VS"""

    P2_distances = []
    for c in P2_list:
        coord = complex_coordinates[c]
        P2_distances.append(distance.euclidean(com_ligand, coord))
    # sort P2list according to distance
    sorted_distance = sorted((x, y) for y, x in enumerate(P2_distances))
    # Choose atom that is closest to Virtual Site
    index_P2 = sorted_distance[0][1]
    P2 = P2_list[index_P2]
    coordinates_P2 = complex_coordinates[P2]

    return coordinates_P2, P2 + 1


def get_virtual_site(coordinates_P1, coordinates_P2, distance_P1_COM):
    """ Get parameters for virtual site"""
    # Vector between two protein atoms
    vec = np.subtract(coordinates_P2, coordinates_P1)
    length = np.math.sqrt(sum(i ** 2 for i in vec))
    # normalize vector
    norm = vec / length
    # Virtual site should be close to COM ligand
    vs = coordinates_P1 + distance_P1_COM * norm
    distance_P1_VS = distance.euclidean(coordinates_P1, vs)
    distance_P2_VS = distance.euclidean(coordinates_P2, vs)
    dist_a = round(distance_P1_VS / (distance_P1_VS + distance_P2_VS), 7)
    dist_b = round(distance_P2_VS / (distance_P1_VS + distance_P2_VS), 7)

    return vs, dist_a, dist_b


def L1_vs(heavy_ligand_coordinates, ligand, vs):
    """ Find ligand atom that is closest to VS and that is restrained to VS with distance restraint"""

    distances = []
    for l in heavy_ligand_coordinates:
        dist = distance.euclidean(l, vs)
        distances.append(round(dist, 3))
    sorted_distance = sorted([x, y] for y, x in enumerate(distances))
    index_ligand = sorted_distance[:1]
    ligand_atom = ligand[index_ligand[0][1]]
    lig_vs = index_ligand[0][0]

    return ligand_atom + 1, lig_vs


def get_contact_restraint_atoms(heavy_ligand_coordinates, protein_coordinates, heavy_ligand, heavy_protein, vs,
                                max_dist=0.35, min_dist_vs=0.4, max_dist_vs=0.7):
    """ Contact restraint between protein and ligand to prevent the ligand from rotating around the center of mass."""
    lig_index = 0
    distances_contre = []
    pl_contre = []
    for c in heavy_ligand_coordinates:
        pro_index = 0
        for p in protein_coordinates:
            d = distance.euclidean(c, p)
            # select a small distance for contact restraint
            if d < max_dist:
                dist_vs = distance.euclidean(c, vs)
                # contact restraint atom should be a little bit away from COM but not too far
                if min_dist_vs < dist_vs < max_dist_vs:
                    distances_contre.append(d)
                    pl_contre.append([heavy_protein[pro_index] + 1, heavy_ligand[lig_index] + 1])
            pro_index += 1
        lig_index += 1
    sorted_distance = sorted([x, y] for y, x in enumerate(distances_contre))

    ###ligand atoms --> +2 (two additional virtual site atoms

    contre_PL = pl_contre[sorted_distance[0][1]]
    dist_contre = round(sorted_distance[0][0], 4)

    return contre_PL, dist_contre


def add_vs_gro(in_gro, out_gro, vs_coordinates):
    """ Add virtual site to gro coordinate file."""

    file = open(in_gro, 'r')
    text = file.readlines()
    file = open(out_gro, 'w')
    co = 0
    # Iterate over complex, insert VS
    prevline = ''
    for idx, line in enumerate(text):
        # Add number of atoms of ligandB to total number of atoms systemA
        if idx == 1:
            line = '%i\n' % (int(line) + 2)
            # Add virtual site
        if 'MOL' in line and co == 0:
            prev = prevline.split()
            vs1 = vs_coordinates[0]
            vs2 = vs_coordinates[1]
            prev_atomnr = int(prev[2]) + 3
            line = '  %iVIR    VS1 %i   %.3f   %.3f   %.3f\n  %iVIR    VS2 %i   %.3f   %.3f   %.3f\n  %s     %s %i   %s   %s   %s\n' \
                   % (int(prev[0][:3]) + 1, int(prev[2]) + 1, vs1[0], vs1[1], vs1[2], \
                      int(prev[0][:3]) + 2, int(prev[2]) + 2, vs2[0], vs2[1], vs2[2], \
                      line.split()[0], line.split()[1], int(prev[2]) + 3, line.split()[3], line.split()[4],
                      line.split()[5])
            print(line)
            co += 1
        if co == 1 and 'VIR' not in line and idx < len(text) - 1:
            atomnr = prev_atomnr + 1
            if atomnr < 10000:
                line = '%s %i   %s' % (line[:15], prev_atomnr + 1, line[23:])
            else:
                line = '%s%i   %s' % (line[:15], prev_atomnr + 1, line[23:])
            prev_atomnr = atomnr
        prevline = line

        file.write(line)
    file.close()
    return


def create_virtual_site(complex_septop, complex_A, lig='MOL'):
    # coordinates septop
    traj = md.load('%s' % complex_septop)
    topology = traj.topology
    # coordinates ligandA: to figure out length ligand A
    traj_1 = md.load('%s' % complex_A)
    topology_1 = traj_1.topology
    ligand_A = topology_1.select('resname %s' % lig).tolist()

    # both ligands
    ligand = topology.select('resname %s' % lig).tolist()
    # Identify ligand A and ligand B
    ligandA = ligand[:len(ligand_A)]
    ligandB = ligand[len(ligand_A):]

    # heavy atoms of protein
    heavy_protein = topology.select('protein and (backbone or name CB)').tolist()
    protein_coordinates = traj.xyz[:, heavy_protein, :]

    vs_coordinates = []
    P1s = []
    P2s = []
    dist_as = []
    ligand_atoms, lig_vs_dist = [], []
    contres_PL, dist_contres = [], []
    # Get one VS for each ligand
    for ligand in [ligandA, ligandB]:
        ligand_traj = traj.atom_slice(ligand, inplace=False)
        # Compute COM of ligand
        com_ligand = md.compute_center_of_mass(ligand_traj)
        topol_ligand = ligand_traj.topology
        heavy_ligand = topol_ligand.select('resname %s and not element H' % lig).tolist()
        heavy_ligand = [ligand[h] for h in heavy_ligand]
        heavy_ligand_coordinates = traj.xyz[:, heavy_ligand, :]

        # Find P1 closest to COM
        coordinates_P1, P1, distance_P1_COM = P1_virtual_site(protein_coordinates[0], heavy_protein, com_ligand[0])
        complex_coordinates = traj.xyz[:, :, :]
        # Find second protein atom that is non colinear with P1 and COM of the ligand (or closest ligand atom)
        P2_list = []
        for r in heavy_protein:
            # restrained_atoms = [P1, L1, r]
            if P1 != r:
                atoms = [complex_coordinates[0][P1], com_ligand[0], complex_coordinates[0][r]]
                # atoms = [complex_coordinates[0][P1], complex_coordinates[0][L1], complex_coordinates[0][r]]
                collinear = _is_collinear(atoms)
                # only choose atoms that are not collinear
                if collinear == True:
                    P2_list.append(r)

        # Take P2 that is closest to COM
        coordinates_P2, P2 = P2_virtual_site(P2_list, complex_coordinates[0], com_ligand[0])

        # Get VS
        vs, dist_a, dist_b = get_virtual_site(coordinates_P1, coordinates_P2, distance_P1_COM)

        # Find heavy ligand atom closest to Virtual site and distance from vs to restrain it to virtual site
        ligand_atom, lig_vs = L1_vs(heavy_ligand_coordinates[0], ligand, vs)

        # Contact restraint to prevent ligand from rotating around virtual site
        contre_PL, dist_contre = get_contact_restraint_atoms(heavy_ligand_coordinates[0], protein_coordinates[0],
                                                             heavy_ligand, heavy_protein, vs, max_dist=0.35,
                                                             min_dist_vs=0.4,
                                                             max_dist_vs=0.7)


        vs_coordinates.append(vs)
        dist_as.append(dist_a)
        P1s.append(P1)
        P2s.append(P2)
        lig_vs_dist.append(lig_vs)
        ligand_atoms.append(ligand_atom)
        contres_PL.append(contre_PL)
        dist_contres.append(dist_contre)

    return vs_coordinates, dist_as, P1s, P2s, lig_vs_dist, ligand_atoms, contres_PL, dist_contres
