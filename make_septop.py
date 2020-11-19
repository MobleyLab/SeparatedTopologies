"""
Functions to build the separated topology of two ligands.
"""

import parmed as pmd

def combine_ligands_top(top_A, top_B, septop, ligand = 'LIG'):
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
    #Choose the protein, water and ions from complex A
    wat = complex_A['HOH', :]
    Cl = complex_A['Cl-', :]
    Na = complex_A['Na+', :]
    prot = complex_A['!(:%s,HOH,Cl-,Na+)'%ligand]

    # Combine different parts
    sep_top = prot + lig1 + lig2 + Na + Cl + wat
    # combine lig1 and lig2 into a single molecule entry
    sep_top.write(septop, [[1, 2]])

    return

def make_section(text):
    #Create dictionary of different section of the topology file

    dic = {}
    start_section = False

    for index, line in enumerate(text):
        # '[ ' marks the start of a section
        if '[ ' in line:
            if not start_section and index==0:
                start_section = True
                key = line
                dic[key] = []
            else:
                #break if already next section
                break
        elif start_section:
            #append all the lines that belong to that section
            dic[key].append(line)

    return dic

def create_A_and_B_state_ligand(line, A_B_state = 'vdwq_q'):
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
    #Turn vdw and electrostatics off
    elif A_B_state == 'vdwq_dummy':
        text = line.split(';')[0] + '   ' + '  d%s  ' % atom_type + '   0.0    ' + mass + '\n'
    #uncharge
    elif A_B_state == 'vdwq_vdw':
        text = line.split(';')[0] + '   ' + '   ' + atom_type + '   0.0    ' + mass + '\n'
    #charge
    elif A_B_state == 'vdw_vdwq':
        text = '   ' + atom_number + '   ' + atom_type + '   ' + residue_nr + '  ' + \
          residue_name + '  ' + atom_name + '  ' + cgnr + '   ' + str(0.0) + '   ' + \
          mass + '    ' + atom_type + '   ' + charge + '   ' + mass + '\n'
    else:
        print('Transformation not implemented yet')

    return text

def create_top(in_top, out_top, A_B_state_ligA, A_B_state_ligB, ligand='LIG'):
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
    while count < end_text:
        #Create dictionary of different sections
        dic = make_section(text[count:])

        count += 1
        # Loop through sections
        for key, value in dic.items():
            #For every atomtype add a dummy-atomtype with no vdW interactions
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
                outtext.append('\n\n')

            #Modify the ligand system
            elif 'atoms' in key and ligand in value[2].split():

                outtext.append(key)
                atomindex_i, atomindex_j = [], []
                for v in value:
                    if v.startswith(';') or v.startswith('\n'):
                        outtext.append(v)
                    #For ligand A create an A and a B state according to A_B_state_ligA
                    if 'LIG' in v and int(v.split()[2]) == 1 and not v.startswith(';'):
                        atomindex_i.append(v.split()[0])
                        line = create_A_and_B_state_ligand(v, A_B_state_ligA)
                        outtext.append(line)
                    #For ligand B create an A and a B state according to A_B_state_ligB
                    if 'LIG' in v and int(v.split()[2]) == 2 and not v.startswith(';'):
                        atomindex_j.append(v.split()[0])
                        line = create_A_and_B_state_ligand(v, A_B_state_ligB)
                        outtext.append(line)
                outtext.append('\n\n[exclusions]\n\n')
                #add exclusions between all atoms of ligand A and ligand B
                for i in atomindex_i:
                    line = '  '.join([i] + atomindex_j)
                    outtext.append('%s\n' % line)
                outtext.append('\n\n')
            else:
                outtext.append(key)
                outtext.append(value)

    for sec in outtext:
        for line in sec:
            file.write(line)
    file.close()

    return out_top