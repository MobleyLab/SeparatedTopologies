"""
Functions to build the separated topology of two ligands.
"""

import parmed as pmd
from SeparatedTopologies import virtual_site as visi


def combine_ligands_top(top_A, top_B, septop, ligand='LIG', water='HOH', Na='Na+', Cl='Cl-'):
    """Create a combined topology file of the complex A and ligand B. Insert ligand B into the same molecule section as ligand A.
    Parameters
    ----------
    top_A : str
        Topology file of complex A
    top_B : str
        Topology file of complex B
    septop : str
        Topology file with ligand B inserted in complex A topology
    """
    complex_A = pmd.load_file(top_A)
    complex_B = pmd.load_file(top_B)

    # Store different molecule_entries separately to be able to assemble the new .top file in the order we want
    lig1 = complex_A[ligand, :]
    lig2 = complex_B[ligand, :]
    # Choose the protein, water and ions from complex A
    wat = complex_A[water, :]
    ClI = complex_A[Cl, :]
    NaI = complex_A[Na, :]
    prot = complex_A['!(:%s,%s,%s,%s)' % (ligand, water, Na, Cl)]

    # Combine different parts
    sep_top = prot + lig1 + lig2 + wat + NaI + ClI
    # combine lig1 and lig2 into a single molecule entry
    sep_top.write(septop, [[1, 2]])

    return


def make_section(text):
    # Create dictionary of different section of the topology file

    dic = {}
    start_section = False

    for index, line in enumerate(text):
        # '[ ' marks the start of a section
        if '[ ' in line:
            if not start_section and index == 0:
                start_section = True
                key = line
                dic[key] = []
            else:
                # break if already next section
                break
        elif start_section:
            # append all the lines that belong to that section
            dic[key].append(line)

    return dic


def create_A_and_B_state_ligand(line, A_B_state='vdwq_q'):
    """Create A and B state topology for a ligand.
    Parameters
    ----------
    line : str
        'Atom line': with atomtype, mass, charge,...
    A_B_state : str
        Interactions in the A state and in the B state.
        vdwq_vdwq: ligand fully interacting in A and B state
        vdwq_vdw: vdw interactions and electrostatics in the A_state, only vdw in the B_state
        vdw_vdwq: charge
        vdw_dummy
        dummy_vdw
        vdwq_dummy
    Returns
    -------
    text : str
        Atoms line for topology file with A and B state parameters
    """
    atom_number = line.split()[0]
    atom_type = line.split()[1]
    residue_nr = line.split()[2]
    residue_name = line.split()[3]
    atom_name = line.split()[4]
    cgnr = line.split()[5]
    charge = line.split()[6]
    mass = line.split()[7]

    # A and B state are the same
    if A_B_state == 'vdwq_vdwq':
        text = line.split(';')[0] + '   ' + atom_type + '   ' + charge + '   ' + mass + '\n'
    # Turn on vdw
    elif A_B_state == 'dummy_vdw':
        charge = str(0.0)
        text = '   ' + atom_number + '   d%s   ' % atom_type + '   ' + residue_nr + '  ' + \
               residue_name + '  ' + atom_name + '   ' + cgnr + '   ' + charge + '   ' + mass + '   ' + \
               atom_type + '   ' + charge + '   ' + mass + '\n'
    # Turn vdw off
    elif A_B_state == 'vdw_dummy':
        charge = str(0.0)
        text = '   ' + atom_number + '   ' + atom_type + '   ' + residue_nr + '   ' + \
               residue_name + '   ' + atom_name + '   ' + cgnr + '   ' + charge + '   ' + mass + \
               '   d%s   ' % atom_type + '   ' + charge + '   ' + mass + '\n'
    # Turn vdw and electrostatics off
    elif A_B_state == 'vdwq_dummy':
        text = line.split(';')[0] + '   ' + '  d%s  ' % atom_type + '   0.0    ' + mass + '\n'
    # uncharge
    elif A_B_state == 'vdwq_vdw':
        text = line.split(';')[0] + '   ' + '   ' + atom_type + '   0.0    ' + mass + '\n'
    # charge
    elif A_B_state == 'vdw_vdwq':
        text = '   ' + atom_number + '   ' + atom_type + '   ' + residue_nr + '  ' + \
               residue_name + '  ' + atom_name + '  ' + cgnr + '   ' + str(0.0) + '   ' + \
               mass + '    ' + atom_type + '   ' + charge + '   ' + mass + '\n'
        # Posre off
    elif A_B_state == 'dummy':
        charge = str(0.0)
        text = '   ' + atom_number + '   d%s   ' % atom_type + '   ' + residue_nr + '  ' + \
               residue_name + '  ' + atom_name + '   ' + cgnr + '   ' + charge + '   ' + mass + '   ' + '\n'
    # Turn vdw and electrostatics off
    elif A_B_state == 'vdwq':
        text = line.split(';')[0] + '\n'
    else:
        print('Transformation not implemented yet')

    return text


def create_top(in_top, out_top, A_B_state_ligA, A_B_state_ligB, complex_septop, complex_A, ligand='LIG'):
    """Create separated topology
    Parameters
    ----------
    in_top : str
        topology file of complex and both ligands (generated from combine_ligands_top)
    out_top: str
        name for output topology file
    A_B_state_ligA : str
        Interactions in the A state and in the B state for ligand A.
        vdwq_vdwq: ligand fully interacting in A and B state
        vdwq_vdw: vdw interactions and electrostatics in the A_state, only vdw in the B_state
        vdw_vdwq: charge
        vdw_dummy
        dummy_vdw
        vdwq_dummy
    A_B_state_ligB: str
        Interactions in the A state and in the B state for ligand B
    ligand = str
        three letter code for the ligand residue, default = 'LIG'
    Examples
    --------
    Turn on vdw ligand B: create_top(top, turnon_vdw_B, 'vdwq_vdwq', 'dummy_vdw')
    Charge ligand B while uncharging ligand A: create_top(top, charge_uncharge, 'vdwq_vdw', 'vdw_vdwq')
    Turn off vdw ligand A: create_top(top, turnoff_vdw_A, 'vdw_dummy', 'vdwq_vdwq')
    """

    file = open(in_top, 'r')
    text = file.readlines()
    file.close()

    file = open(out_top, 'w')
    end_text = len(text)
    count = 0
    outtext = []
    section = 0
    while count < end_text:
        # Create dictionary of different sections
        dic = make_section(text[count:])

        count += 1

        # Loop through sections
        for key, value in dic.items():
            # For every atomtype add a dummy-atomtype with no vdW interactions
            if 'atomtypes' in key:
                outtext.append(key)
                outtext.append(value)
                for v in value:
                    if v.startswith(';') or v.startswith('\n'):
                        continue
                    else:
                        at_line = v.split()
                        dummy = 'd' + at_line[0] + '       ' + at_line[1] + '     ' + at_line[
                            2] + ' 0.0 ' + ' A ' + ' 0.0 ' + ' 0.0\n'
                        outtext.append(dummy)
                outtext.append('VS    VS  0.0000 0.0000  D  0.0 0.0')
                outtext.append('\n\n')


            # Modify protein system to include VS
            elif 'atoms' in key and section == 3:
                outtext.append(key)
                outtext.append(value)
                prev = value[-2].split()
                vs1 = int(prev[0]) + 1
                vs2 = int(prev[0]) + 2
                outtext.append(' %i  VS  %i  VIR   VS1   %i  0.000 0.000 \n %i  VS  %i  VIR   VS2   %i  0.000 0.000 \n' \
                               % (vs1, int(prev[2]) + 1, vs1, vs2, int(prev[2]) + 2, vs2))

            elif 'dihedral' in key and section == 7:
                outtext.append(key)
                outtext.append(value)
                vs_coord, dist_as, P1s, P2s, lig_vs_dist, ligand_atoms, contres_PL, dist_contres = visi.create_virtual_site(
                    complex_septop, complex_A, lig=ligand)
                outtext.append('[ virtual_sites2 ]\n %i   %i %i 1 %.7f\n %i   %i %i 1 %.7f\n\n' \
                               % (vs1, P1s[0], P2s[0], dist_as[0], vs2, P1s[1], P2s[1], dist_as[1]))
            # Modify the ligand system
            elif 'atoms' in key and ligand in value[2].split():

                outtext.append(key)
                atomindex_i, atomindex_j = [], []
                for v in value:
                    if v.startswith(';') or v.startswith('\n'):
                        outtext.append(v)
                    # For ligand A create an A and a B state according to A_B_state_ligA
                    if ligand in v and int(v.split()[2]) == 1 and not v.startswith(';'):
                        atomindex_i.append(v.split()[0])
                        line = create_A_and_B_state_ligand(v, A_B_state_ligA)
                        outtext.append(line)
                    # For ligand B create an A and a B state according to A_B_state_ligB
                    if ligand in v and int(v.split()[2]) == 2 and not v.startswith(';'):
                        atomindex_j.append(v.split()[0])
                        line = create_A_and_B_state_ligand(v, A_B_state_ligB)
                        outtext.append(line)
                outtext.append('\n\n[exclusions]\n\n')
                # add exclusions between all atoms of ligand A and ligand B
                for i in atomindex_i:
                    line = '  '.join([i] + atomindex_j)
                    outtext.append('%s\n' % line)
                outtext.append('\n\n')
            else:
                outtext.append(key)
                outtext.append(value)

            section += 1

    for sec in outtext:
        for line in sec:
            file.write(line)

    file.write('\n[ intermolecular_interactions ]\n[ bonds ]\n; ai    aj    type   bA    kA    bB    kB\n')
    fc_A = 8368.0
    fc_B = 8368.0
    if A_B_state_ligA == 'dummy':
        fc_B = 0.0
    if A_B_state_ligB == 'dummy':
        fc_A = 0.0

    file.write('%i   %i   6   %.3f    %.1f   %.3f   %.1f\n' % (
    vs1, ligand_atoms[0] + 2, lig_vs_dist[0], fc_A, lig_vs_dist[0], fc_B))
    file.write('%i   %i   6   %.3f    %.1f   %.3f   %.1f\n' % (
    vs2, ligand_atoms[1] + 2, lig_vs_dist[1], fc_A, lig_vs_dist[1], fc_B))
    file.write('%i   %i   6   %.3f    %.1f   %.3f   %.1f\n' % (
    contres_PL[0][0], contres_PL[0][1] + 2, dist_contres[0], fc_A, dist_contres[0], fc_B))
    file.write('%i   %i   6   %.3f    %.1f   %.3f   %.1f\n' % (
    contres_PL[1][0], contres_PL[1][1] + 2, dist_contres[1], fc_A, dist_contres[1], fc_B))

    file.close()

    return out_top, vs_coord
