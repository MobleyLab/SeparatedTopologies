import parmed as pmd
from SeparatedTopologies import solvent as so

#Name of project directory
dir = ''
#List of ligands
compounds = []
#Three letter code ligands
lig = 'UNL'

for f in compounds:
    path = '%s/%s' % (dir,f)

    top = '%s/solvent.top'%path
    so.create_top(top, top, 'vdwq_dummy', ligand=lig)
