"""
Functions to build the separated topology of two ligands.
"""

import parmed as pmd

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


def create_A_and_B_state_ligand(line, A_B_state='vdwq_q', lig = 1):
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
        vdwq_scaled-vdw
        scaled-vdw_dummy
        dummy_scaled-vdwq
        scaled-vdwq_vdwq
    lig: int
        if first or second ligand
    Returns
    -------
    text : str
        Atoms line for topology file with A and B state parameters
    """
    atom_number = line.split()[0]
    atom_type = line.split()[1]
    if lig == 1:
        atom_type = "".join(('LIG1_', atom_type))
    if lig == 2:
        atom_type = "".join(('LIG2_', atom_type))
    scaled_atomtype = "".join(('scaled_', atom_type))
    residue_nr = line.split()[2]
    residue_name = line.split()[3]
    atom_name = line.split()[4]
    cgnr = line.split()[5]
    charge = line.split()[6]
    mass = line.split()[7]

    # A and B state are the same
    if A_B_state == 'vdwq_vdwq':
        text = '   ' + atom_number + '  ' + atom_type + '   ' + residue_nr + '  ' + \
               residue_name + '  ' + atom_name + '   ' + cgnr + '   ' + charge + '   ' + mass + '   ' + \
               atom_type + '   ' + charge + '   ' + mass + '\n'
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
        text = '   ' + atom_number + '   ' + atom_type + '   ' + residue_nr + '   ' + \
               residue_name + '   ' + atom_name + '   ' + cgnr + '   ' + charge + '   ' + mass + \
               '   d%s   ' % atom_type + '   ' + '   0.0    ' + '   ' + mass + '\n'
    # uncharge
    elif A_B_state == 'vdwq_vdw':
        text = '   ' + atom_number + '  ' + atom_type + '   ' + residue_nr + '  ' + \
               residue_name + '  ' + atom_name + '   ' + cgnr + '   ' + charge + '   ' + mass + '   ' + \
               atom_type + '   ' + '   0.0    ' + '   ' + mass + '\n'
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
        text = '   ' + atom_number + '  ' + atom_type + '   ' + residue_nr + '  ' + \
               residue_name + '  ' + atom_name + '   ' + cgnr + '   ' + charge + '   ' + mass + '   ' + '\n'
    # vdw till gamma and charges
    elif A_B_state == 'vdwq_scaled-vdw':
        text = '   ' + atom_number + '  ' + atom_type + '   ' + residue_nr + '  ' + \
               residue_name + '  ' + atom_name + '   ' + cgnr + '   ' + charge + '   ' + mass + '   ' + \
               scaled_atomtype + '   ' + '   0.0    ' + '   ' + mass + '\n'
        # Turn vdw off from scaled ligand
    elif A_B_state == 'scaled-vdw_dummy':
        charge = str(0.0)
        text = '   ' + atom_number + '   ' + scaled_atomtype + '   ' + residue_nr + '   ' + \
               residue_name + '   ' + atom_name + '   ' + cgnr + '   ' + charge + '   ' + mass + \
               '   d%s   ' % atom_type + '   ' + charge + '   ' + mass + '\n'
        # Turn on vdw
    elif A_B_state == 'dummy_scaled-vdwq':
        text = '   ' + atom_number + '   d%s   ' % atom_type + '   ' + residue_nr + '  ' + \
               residue_name + '  ' + atom_name + '   ' + cgnr + '   ' + str(0.0) + '   ' + mass + '   ' + \
               scaled_atomtype + '   ' + charge + '   ' + mass + '\n'
        # turn on rest of vdw
    elif A_B_state == 'scaled-vdwq_vdwq':
        text = '   ' + atom_number + '   ' + scaled_atomtype + '   ' + residue_nr + '  ' + \
                residue_name + '  ' + atom_name + '   ' + cgnr + '   ' + charge + '   ' + mass + '   ' + \
                atom_type + '   ' + charge + '   ' + mass + '\n'
    else:
        print('Transformation not implemented yet')

    return text

def atom_types_ligand(in_top, ligand='LIG'):
    """Store atom types ligand in a list. Atom types have to come from individual ligand files, because parmed just
    combines different atom types to a single one if they have the same name.
    Parameters
    ----------
    in_top : str
        topology file of complex and both ligands (generated from combine_ligands_top)
    ligand = str
        three letter code for the ligand residue, default = 'LIG'
    """

    file = open(in_top, 'r')
    text = file.readlines()
    file.close()
    end_text = len(text)
    # at = []
    atomtype = []
    for line in text:

        if ligand in line and not line.startswith(';') and line.split()[0].isdigit():
            at = line.split()[1]
            for l in text:
                if l.startswith(';') or l.startswith('\n'):
                    continue
                if l.split()[0] == at and l not in atomtype:
                    atomtype.append(l)

    return atomtype

def create_top(in_top, out_top, gamma, A_B_state_ligA, A_B_state_ligB, in_top_A, in_top_B, ligand='LIG'):
    """Create separated topology
    Parameters
    ----------
    in_top : str
        topology file of complex and both ligands (generated from combine_ligands_top)
    out_top: str
        name for output topology file
    gamma: int
        scaling parameter for REST scaling
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

    atomtype_i = atom_types_ligand(in_top_A, ligand=ligand)
    atomtype_j = atom_types_ligand(in_top_B, ligand=ligand)
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
                for v in value:
                    if v.startswith(';') or v.startswith('\n'):
                        continue
                    if v.split()[0] in [l.split()[0] for l in atomtype_i]:
                        continue
                    if v.split()[0] in [l.split()[0] for l in atomtype_j]:
                        continue
                    else:
                        outtext.append(v)
                for atomtype in [atomtype_i, atomtype_j]:
                    for l in atomtype:
                        if atomtype == atomtype_i:
                            new_at = "".join(('LIG1_', l))
                        elif atomtype == atomtype_j:
                            new_at = "".join(('LIG2_', l))
                        at_line = new_at.split()
                        #Scale the epsilon of the vdw terms by a constant gamma
                        scaled_at = 'scaled_' + at_line[0] + '       ' + at_line[1] + '     ' + at_line[
                                2] + ' 0.0 ' + ' A ' + at_line[5] + '  ' + str(float(at_line[6])*gamma) + '\n'

                        ###To do: decide if I only allow those 4 endstates or make it more flexible!!!

                        dummy = 'd' + at_line[0] + '       ' + at_line[1] + '     ' + at_line[
                                2] + ' 0.0 ' + ' A ' + ' 0.0 ' + ' 0.0\n'

                        outtext.append(new_at)
                        outtext.append(scaled_at)
                        outtext.append(dummy)

                outtext.append('\n\n')

                outtext.append('[ nonbond_params ]\n')
                for i in atomtype_i:
                    i = "".join(('LIG1_', i.split()[0]))
                    scaled_i = "".join(('scaled_', i.split()[0]))
                    for j in atomtype_j:
                        j = "".join(('LIG2_', j.split()[0]))
                        scaled_j = "".join(('scaled_', j.split()[0]))
                        nb = i + '  ' + j + '   1   0   0\n'
                        scaled_nb = scaled_i + '  ' + scaled_j + '   1   0   0\n'
                        outtext.append(nb)
                        outtext.append(scaled_nb)

                outtext.append('\n\n')


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
                        line = create_A_and_B_state_ligand(v, A_B_state_ligA, lig = 1)
                        outtext.append(line)
                    # For ligand B create an A and a B state according to A_B_state_ligB
                    if ligand in v and int(v.split()[2]) == 2 and not v.startswith(';'):
                        atomindex_j.append(v.split()[0])
                        line = create_A_and_B_state_ligand(v, A_B_state_ligB, lig = 2)
                        outtext.append(line)

                outtext.append('\n\n')
            else:
                outtext.append(key)
                outtext.append(value)

            section += 1

    for sec in outtext:
        for line in sec:
            file.write(line)

    file.close()

    return out_top
