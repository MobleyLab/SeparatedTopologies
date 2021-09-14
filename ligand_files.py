from openeye import oechem

def ligand_sdf(pdb, outfile):
    ''' Get ligand sdf from complex pdb '''

    complex_A = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    ofs = oechem.oemolostream()
    ofs.open(outfile)
    ifs.SetFlavor(oechem.OEFormat_PDB,
                  oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)
    ifs.open(pdb)
    oechem.OEReadMolecule(ifs, complex_A)
    ifs.close()

    lig = oechem.OEGraphMol()
    prot = oechem.OEGraphMol()
    wat = oechem.OEGraphMol()
    other = oechem.OEGraphMol()
    #get ligand
    oechem.OESplitMolComplex(lig, prot, wat, other, complex_A)

    oechem.OEWriteMolecule(ofs, lig)

    return