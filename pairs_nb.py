import math

compound_A = 'lig_ejm_55'
compound_B = 'lig_ejm_54'
edge_A_B = '../edge_ejm_55_ejm_54'
lig = 'MOL'
top_A = '%s/complex/complex.top'%compound_A
top_B = '%s/complex/complex.top'%compound_B
top = '%s/complex.top'%edge_A_B

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

def create_top(in_top, A_B_state_ligA, A_B_state_ligB, ligand='LIG'):
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

    end_text = len(text)
    count = 0
    outtext = []
    section = 0
    fake_bond = []
    while count < end_text:
        # Create dictionary of different sections
        dic = make_section(text[count:])

        count += 1
        # Loop through sections
        for key, value in dic.items():
            # Modify the ligand system
            if 'atoms' in key and ligand in value[2].split():
                atomindex_i, atomindex_j = [], []
                charges_i, charges_j = [], []
                for v in value:

                    # For ligand A create an A and a B state according to A_B_state_ligA
                    if ligand in v and int(v.split()[2]) == 1 and not v.startswith(';'):
                        atomindex_i.append(v.split()[0])
                        charges_i.append(v.split()[6])
                    # For ligand B create an A and a B state according to A_B_state_ligB
                    if ligand in v and int(v.split()[2]) == 2 and not v.startswith(';'):
                        atomindex_j.append(v.split()[0])
                        charges_j.append(v.split()[6])
                # add exclusions between all atoms of ligand A and ligand B
                print(atomindex_i)
                print(atomindex_j)
                for ind_i, i in enumerate(atomindex_i):
                    for ind_j, j in enumerate(atomindex_j):
                        line = '%i   %i   1   %.8f    %.8f    0.0    0.0'%(int(i), int(j), -(float(charges_i[ind_i])), -(float(charges_j[ind_j])))
                        fake_bond.append('%s\n' % line)
    file = open('pairs_nb.txt', 'w')
    for sec in fake_bond:
        for line in sec:
            file.write(line)

    file.close()
    return

def create_top_charge(in_top, lam, A_B_state_ligA, A_B_state_ligB, outfile, ligand='LIG'):
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

    end_text = len(text)
    count = 0
    outtext = []
    section = 0
    fake_bond = []
    while count < end_text:
        # Create dictionary of different sections
        dic = make_section(text[count:])

        count += 1
        # Loop through sections
        for key, value in dic.items():
            # Modify the ligand system
            if 'atoms' in key and ligand in value[2].split():
                atomindex_i, atomindex_j = [], []
                charges_i, charges_j = [], []
                for v in value:

                    # For ligand A create an A and a B state according to A_B_state_ligA
                    if ligand in v and int(v.split()[2]) == 1 and not v.startswith(';'):
                        atomindex_i.append(v.split()[0])
                        charges_i.append(v.split()[6])
                    # For ligand B create an A and a B state according to A_B_state_ligB
                    if ligand in v and int(v.split()[2]) == 2 and not v.startswith(';'):
                        atomindex_j.append(v.split()[0])
                        charges_j.append(v.split()[6])
                # add exclusions between all atoms of ligand A and ligand B
                c_i = [float(i) for i in charges_i]
                c_j = [float(i) for i in charges_j]
                cn_i = [float(i)*math.sqrt(1-float(lam)) for i in charges_i]
                cn_j = [float(i)*math.sqrt(float(lam)) for i in charges_j]
                print(sum(c_i))
                print(sum(c_j))
                print(sum(cn_i))
                print(sum(cn_j))
                print(atomindex_i)
                print(atomindex_j)
                for ind_i, i in enumerate(atomindex_i):
                    for ind_j, j in enumerate(atomindex_j):
                        line = '%i   %i   1   %s    %s    0.0    0.0'%(int(i), int(j), str(round(-float(charges_i[ind_i])*math.sqrt(1-float(lam)),8)), str(round(float(charges_j[ind_j])*math.sqrt(float(lam)),8)))
                        fake_bond.append('%s\n' % line)
    file = open(outfile, 'w')
    for sec in fake_bond:
        for line in sec:
            file.write(line)

    file.close()
    return

edge_A_B = 'edge_ejm_55_ejm_54_solvent'
#lambdas = ['0.0','0.25', '0.5', '0.75','1.0']
lambdas = ['0.25', '0.5', '0.75']
for ind,l in enumerate(lambdas):
    create_top_charge(top, l, 'vdwq', 'vdwq', '../%s/pairs_nb_%i.itp'%(edge_A_B,ind+1),ligand='MOL')