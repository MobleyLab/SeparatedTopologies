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

def ligand_pdb(sdf, outfile):
    ''' Get ligand sdf from complex pdb '''

    complex_A = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    ifs.SetFormat(oechem.OEFormat_SDF)
    ofs = oechem.oemolostream()
    ofs.open(outfile)
    # flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_TER
    # ofs.SetFlavor(oechem.OEFormat_PDB, flavor)
    ifs.open(sdf)
    oechem.OEReadMolecule(ifs, complex_A)
    ifs.close()

    # oechem.OEWriteMolecule(ofs, complex_A)
    oechem.OEWritePDBFile(ofs, complex_A)

    return

def ligand_mol2(sdf, outfile):
    ''' Get ligand sdf from complex pdb '''

    complex_A = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    ofs = oechem.oemolostream()
    ofs.open(outfile)
    # ofs.SetFlavor(oechem.OEFormat_PDB,
    #               oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)
    ifs.open(sdf)
    oechem.OEReadMolecule(ifs, complex_A)
    ifs.close()

    oechem.OEWriteMolecule(ofs, complex_A)

    return