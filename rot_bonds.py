from openeye import oechem
import mdtraj as md
import numpy as np

def rota_bonds(file):


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

    if oechem.OESplitMolComplex(lig, prot, wat, other, complex_A):
        rot_bonds = []
        for atom in lig.GetAtoms():

            for bond in atom.GetBonds():
                rot_bond = []
                if bond.IsRotor():
                    nbor = bond.GetNbr(atom)
                    if nbor.GetIdx() > atom.GetIdx():
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


def get_dihedrals(file, lig='MOL'):
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
        dih.append([r + 3 for r in rb])
        values.append(round(float(dih1[0]), 2))
    return dih, values, len_lig

fileA = 'lig_ejm_44/complex/ions.pdb'
fileB = 'lig_ejm_42/complex/ions.pdb'

dih, values, len_ligA = get_dihedrals(fileA)
print(dih)
print(values)
dih, values, len_lig = get_dihedrals(fileB)
dih_B = []
for d in dih:
    d = [x+len_ligA for x in d]
    dih_B.append(d)

print(dih_B)

