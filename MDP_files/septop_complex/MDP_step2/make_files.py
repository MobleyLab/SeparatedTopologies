#!/bin/env python

import shutil
import os

outdir = 'complete_files'
mdpfiles = ['em.X.mdp', 'nvt.X.mdp', 'prod.X.mdp']
lambdanr=8

if not os.path.isdir(outdir): os.mkdir(outdir)
for i in range(0, lambdanr):
    # Copy and edit mdp files
    for mdpfile in mdpfiles:
        mdpfilename = os.path.join( outdir, mdpfile.replace('.X.', '.%s.' % i))
        shutil.copy(mdpfile, mdpfilename)
        file = open(mdpfilename, 'r')
        text = file.readlines()
        file.close()
        file = open(mdpfilename, 'w')
        for idx, line in enumerate(text):
           if 'init_lambda_state' in line and line[0] != ';':
                 line = 'init_lambda_state        = %s\n' % i
           file.write(line)
        file.close()

