from openeye import oechem, oespruce
import parmed as pmd
import shutil
import mdtraj as md

def align_complexes(pdb_A, pdb_B, out_B):
    """Align 2 structures with oespruce, Structure A is the reference
    Parameters
    ----------
    pdb_A : str
       pdb structure Complex A
    pdb_B : str
        pdb structure Complex B
    out_B : str
        pdb structure, Aligned structure of complex B
    """
    complex_A = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB,
                  oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)
    ifs.open(pdb_A)
    oechem.OEReadMolecule(ifs, complex_A)
    ifs.close()

    complex_B = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB,
                  oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)
    ifs.open(pdb_B)
    oechem.OEReadMolecule(ifs, complex_B)
    ifs.close()

    # superposition = oespruce.OEStructuralSuperposition(ref_prot, fit_prot)
    superposition = oespruce.OEStructuralSuperposition(complex_A, complex_B)
    # superposition.Transform(fit_prot)
    superposition.Transform(complex_B)

    ofs = oechem.oemolostream()
    flavor = ofs.GetFlavor(oechem.OEFormat_PDB) ^ oechem.OEOFlavor_PDB_OrderAtoms

    ofs.SetFlavor(oechem.OEFormat_PDB, flavor)

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
    print(pdb)
    pdb.save(gro, overwrite=True)

def gro2pdb(gro, pdb):
    #Takes a pdb file and converts it to a .gro file
    gro = pmd.load_file(gro)
    gro.save(pdb, overwrite=True)

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

def edit_indices(in_gro, out_gro):
    """ Edit indices gro coordinate file."""

    file = open(in_gro, 'r')
    text = file.readlines()
    file = open(out_gro, 'w')
    co = 0
    # Iterate over complex
    prevline = ''
    count = 0
    for idx, line in enumerate(text):
        if 'MOL' in line and count == 0:
            count = 1
            first_MOL_line = idx
            prev = prevline.split()
            prev_atomnr = int(prev[2])
            atomnr = prev_atomnr + 1
            if atomnr < 10000:
                line = '%s %i   %s' % (line[:15], prev_atomnr + 1, line[23:])
            else:
                line = '%s%i   %s' % (line[:15], prev_atomnr + 1, line[23:])
        #if count >= 1 and atomnr != prev_atomnr + 1 and idx < len(text) - 1:
        if count >= 1 and idx > first_MOL_line and idx < len(text) - 1:
            atomnr = atomnr + 1
            if atomnr < 10000:
                line = '%s %i   %s' % (line[:15], atomnr, line[23:])
            else:
                line = '%s%i   %s' % (line[:15], atomnr, line[23:])
        prevline = line

        file.write(line)
    file.close()
    return

def edit_indices_solvent(in_gro, out_gro):
    """ Edit indices gro coordinate file."""

    file = open(in_gro, 'r')
    text = file.readlines()
    file = open(out_gro, 'w')
    co = 0
    # Iterate over complex
    prevline = ''
    count = 0
    for idx, line in enumerate(text):
        if 'MOL' in line and count == 0:
            count = 1
            first_MOL_line = idx
            prev = line.split()
            prev_atomnr = int(prev[2])
            atomnr = prev_atomnr
            if atomnr < 10000:
                line = '%s %i   %s' % (line[:15], prev_atomnr + 1, line[23:])
            else:
                line = '%s%i   %s' % (line[:15], prev_atomnr + 1, line[23:])
        #if count >= 1 and atomnr != prev_atomnr + 1 and idx < len(text) - 1:
        if count >= 1 and idx > first_MOL_line and idx < len(text) - 1:
            atomnr = atomnr + 1
            if atomnr < 10000:
                line = '%s %i   %s' % (line[:15], atomnr, line[23:])
            else:
                line = '%s%i   %s' % (line[:15], atomnr, line[23:])
        prevline = line

        file.write(line)
    file.close()
    return