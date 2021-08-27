from openeye import oechem
import mdtraj as md
import numpy as np
from openeye.oechem import *
from openeye.oedepict import *

def rota_bonds(file):
    ''' Find rotatable bonds ligand. Check for symmetry and only include non-symmetric rotatable bonds '''

    complex_A = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB,
                  oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)
    ifs.open(file)
    oechem.OEReadMolecule(ifs, complex_A)
    ifs.close()

    lig = oechem.OEGraphMol()
    prot = oechem.OEGraphMol()
    wat = oechem.OEGraphMol()
    other = oechem.OEGraphMol()
    #get ligand
    if oechem.OESplitMolComplex(lig, prot, wat, other, complex_A):
        rot_bonds = []
        OEPerceiveSymmetry(lig)

        for atom in lig.GetAtoms():
            for bond in atom.GetBonds():
                rot_bond = []
                #check if bond is rotatable
                if bond.IsRotor():
                    #find neighboring atom (atom2)
                    nbor = bond.GetNbr(atom)
                    #only count bond once (we also loop over neighboring atom, same bond)
                    if nbor.GetIdx() > atom.GetIdx():
                        symmetry = False
                        # loop over all atoms that atom1 is connected with, get their symmetry class
                        sym_a1 = []
                        for a in atom.GetAtoms():
                            sym = a.GetSymmetryClass()
                            #if the same symmetry group appears more than once: there is symmetry
                            if sym in sym_a1:
                                symmetry = True
                            sym_a1.append(sym)
                        # loop over all atoms that atom2 is connected with, get their symmetry class
                        # do this since symmetry can occur at either side of the rotatable bond
                        sym_a2 = []
                        for a in nbor.GetAtoms():
                            sym = a.GetSymmetryClass()
                            if sym in sym_a2:
                                symmetry = True
                            sym_a2.append(sym)
                        # if the rotatable bond is not symmetric: save the 4 atoms of that dihedral
                        if symmetry == False:
                            for a in atom.GetAtoms():
                                if a.GetIdx() != nbor.GetIdx() and a.IsHydrogen() == False:
                                    d1 = a.GetIdx()

                            for n in nbor.GetAtoms():
                                if n.GetIdx() != atom.GetIdx() and n.IsHydrogen() == False:
                                    d4 = n.GetIdx()
                            rot_bond.append(d1)
                            rot_bond.append(atom.GetIdx())
                            rot_bond.append(nbor.GetIdx())
                            rot_bond.append(d4)
                            rot_bonds.append(rot_bond)

    return rot_bonds


def get_dihedrals(file, lig):
    '''Get dihedral around rotatable bond'''
    traj = md.load(file)
    topology = traj.topology
    ligand = topology.select('resname %s' % lig).tolist()
    len_lig = len(ligand)
    rot_bonds = rota_bonds(file)
    dih = []
    values = []
    for rb in rot_bonds:
        rb = [ligand[r] for r in rb]
        dih1 = np.rad2deg(md.compute_dihedrals(traj, [np.array(rb)]))
        dih.append([r + 1 for r in rb])
        values.append(round(float(dih1[0]), 2))
    return dih, values, len_lig

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

def restrain_rot_bonds(fileA, fileB, folder,lig = 'MOL'):
    dih_A, values_A, len_ligA = get_dihedrals(fileA, lig)
    dih, values_B, len_lig = get_dihedrals(fileB, lig)
    dih_B = []
    for d in dih:
        d = [x + len_ligA for x in d]
        dih_B.append(d)

    # dih = dih_A + dih_B
    # values = values_A + values_B
    write_itp_restraints(dih_A, values_A, 5, 5, '%s/rot_bonds_A.itp'%folder)
    write_itp_restraints(dih_B, values_B, 5, 5, '%s/rot_bonds_B.itp' % folder)
    write_itp_restraints(dih_A, values_A, 0, 5, '%s/rot_bonds_A_on.itp' % folder)
    write_itp_restraints(dih_B, values_B, 5, 0, '%s/rot_bonds_B_off.itp' % folder)

    return

def include_itp_in_top(top, idpfile):
    with open(top) as file:
        newline = ''
        for line in file:
            for part in line.split():
                if 'rot_bonds' not in part:
                    newline = '\n#include "%s"'%idpfile
        file.close()
    file = open(top, 'a')
    file.write(newline)
    file.close()