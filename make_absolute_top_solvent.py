import parmed as pmd
import ligand_files as lf
import rot_bonds as rotbond
import mdtraj as md
import numpy as np

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

def create_top(in_top, out_top, A_B_state_ligA, ligand='LIG'):
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
                        line = create_A_and_B_state_ligand(v, A_B_state_ligA)
                        outtext.append(line)
            else:
                outtext.append(key)
                outtext.append(value)

            section += 1

    for sec in outtext:
        for line in sec:
            file.write(line)
    file.close()

    return out_top

def write_itp_restraints(dih, values, forceconst_A, forceconst_B, file):
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

    fc_rad_a = forceconst_A * 4.184
    fc_rad_b = forceconst_B * 4.184

    file = open(file, 'w')
    file.write('[ intermolecular_interactions ] \n[ dihedrals ] \n')
    file.write('; ai     aj    ak    al    type     thA      fcA       thB      fcB\n')
    for inx, d in enumerate(dih):
        file.write(' %s   %s   %s   %s   2   %.2f   %.2f   %.2f   %.2f\n' % (
        d[0], d[1], d[2], d[3], values[inx], fc_rad_a, values[inx], fc_rad_b))

    file.close()

    return

def get_dihedrals(ligand,solvent,complex, lig):
    '''Get dihedral around rotatable bond'''
    traj = md.load(complex)
    topology = traj.topology
    ligand_top = topology.select('resname %s' % lig).tolist()
    len_lig = len(ligand_top)
    print(len_lig)
    rot_bonds = rotbond.rota_bonds(ligand)
    traj_solvent = md.load(solvent)
    topology_solvent = traj_solvent.topology
    ligand_top_solvent = topology_solvent.select('resname UNL').tolist()
    print(len(ligand_top_solvent))
    dih = []
    values = []
    for rb in rot_bonds:
        rb_solvent = [ligand_top_solvent[r] for r in rb]
        rb = [ligand_top[r] for r in rb]
        dih1 = np.rad2deg(md.compute_dihedrals(traj, [np.array(rb)]))
        dih.append([r + 1 for r in rb_solvent])
        values.append(round(float(dih1[0]), 2))
    return dih, values, len_lig

def restrain_rot_bonds(ligand, solvent,pdb,lig,folder):
    dih_A, values_A, len_ligA = get_dihedrals(ligand,solvent,pdb,lig)
    write_itp_restraints(dih_A, values_A, 0, 5, '%s/rot_bonds_on.itp' % folder)

    return

