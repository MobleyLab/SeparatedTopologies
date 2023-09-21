import parmed as pmd
import make_absolute_top_solvent as so

#Name of project directory
dir = ''
#List of ligands
compounds = []
#Three letter code ligands
lig = 'UNL'

for f in compounds:
    path = '%s/%s' % (dir,f)
    #specify output file name
    top = '%s/solvent.top'%path
    #Add dummy B state to the topology file
    so.create_top(top, top, 'vdwq_dummy', ligand=lig)

