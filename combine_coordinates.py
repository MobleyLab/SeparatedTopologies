from openeye import oechem, oespruce
import parmed as pmd
import shutil
import mdtraj as md

def align_complexes_mdtraj(pdb_A, pdb_B, out_B):
    """Align 2 structures with oespruce, Structure A is the reference
    Parameters
    ----------
    pdb_A : str
        pdb structure Complex A
    pdb_B : str
        pdb structure Complex B
    out_B : str
        pdb structure Aligned structure of complex B
    """

    ligand_A = md.load(pdb_A)
    ligand_B = md.load(pdb_B)
    ligand_B.superpose(ligand_A)
    ligand_B.save(out_B)

    return

def align_complexes(pdb_A, pdb_B, out_B):
    """Align 2 structures with oespruce, Structure A is the reference
    Parameters
    ----------
    pdb_A : str
        pdb structure Complex A
    pdb_B : str
        pdb structure Complex B
    out_B : str
        pdb structure Aligned structure of complex B
    """

    complex_A = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB,
                  oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)
    ifs.open(pdb_A)
    oechem.OEReadMolecule(ifs, complex_A)
    # oechem.OEReadPDBFile(ifs, complex_A)
    ifs.close()

    complex_B = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB,
                  oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)
    ifs.open(pdb_B)
    oechem.OEReadMolecule(ifs, complex_B)
    # oechem.OEReadPDBFile(ifs, complex_B)
    ifs.close()

    # superposition = oespruce.OEStructuralSuperposition(ref_prot, fit_prot)
    superposition = oespruce.OEStructuralSuperposition(complex_A, complex_B)

    # superposition.Transform(fit_prot)
    superposition.Transform(complex_B)

    ofs = oechem.oemolostream()
    ofs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)
    ofs.open(out_B)
    oechem.OEWriteMolecule(ofs, complex_B)
    ofs.close()

    return


def combine_ligands_gro(in_file_A, in_file_B, out_file, ligand_A='MOL', ligand_B='MOL'):
    """Add ligand B coordinates to coordinate (.gro) file of ligand A in complex with protein
    Parameters
    ----------
    in_file_A : .gro file
        Complex A coordinate file
    in_file_B : .gro file
        Complex B coordinate file
    out_file : .gro file
        Complex A + ligand B coordinate file
    ligand_A: str
        Three letter code for ligand A
    ligand_B: str
        Three letter code for ligand B
    """
    # make copy of in_file_A ---> new outfile
    shutil.copy(in_file_A, out_file)
    file = open(in_file_B, 'r')
    text_B = file.readlines()
    file.close()
    file = open(out_file, 'r')
    text_A = file.readlines()
    file = open(out_file, 'w')
    lig2 = []
    #Iterate over complex B, store ligand B lines
    for idx, line in enumerate(text_B):
        if ligand_B in line:
            lig2.append(line)
    print(lig2)

    lig1 = []
    # Iterate over complex A, store ligand A lines
    for idx, line in enumerate(text_A):
        if ligand_A in line:
            lig1.append(line)
    # Iterate over complex A, insert ligand B lines
    for idx, line in enumerate(text_A):
        # Add number of atoms of ligandB to total number of atoms systemA
        if idx == 1:
            line = '%i\n' % (int(line) + len(lig2))
            # Insert ligandB coordinates right after ligandA
        if line == lig1[-1]:
            count = 0
            for i in lig2:
                index = idx + count + 1
                text_A.insert(index, i)
                count += 1
        file.write(line)
    file.close()

    return out_file

def pdb2gro(pdb, gro):
    #Takes a pdb file and converts it to a .gro file
    pdb = pmd.load_file(pdb)
    pdb.save(gro, overwrite=True)

def ligand_heavyatoms_ndx(traj, ligand='LIG'):
    """Write index file with ligand heavy atoms for position restraints.
    Parameters
    ----------
    traj : mdtraj trajectory
        Mdtraj object with coordinates of the system (e.g. from .gro file)
    ligand : str
        Three letter code for ligand name
    Returns
    -------

    """
    traj = md.load(traj)
    topology = traj.topology
    ligand = topology.select('resname %s' % ligand).tolist()
    ligand_traj = traj.atom_slice(ligand, inplace=False)
    topology_ligand = ligand_traj.topology
    heavy_ligand = topology_ligand.select('resname LIG and not element H').tolist()
    heavy_ligand = [i+1 for i in heavy_ligand]

    return heavy_ligand

# pdb_A = 'cpd1/complex.pdb'
# pdb_B = 'cpd6/complex.pdb'
#
# fit_B = 'cpd6/complex_fit.pdb'
# gro_B = 'cpd6/complex_fit.gro'
# gro_A = 'cpd1/complex.gro'
# complex = 'complex.gro'
# align_complexes(pdb_A, pdb_B, fit_B)
# #Convert .pdb to .gro
# pdb2gro(fit_B, gro_B)
#
# combine_ligands_gro(gro_A, gro_B, complex)