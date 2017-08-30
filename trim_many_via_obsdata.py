#!/usr/bin/env python
"""
Code to create many trimmed model grids for batch runs
  Saves time by only reading the huge modelsed grid once
"""

# system imports
from __future__ import print_function
import sys
import os
import argparse
import time

# BEAST imports
import datamodel_big as datamodel
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel 
from beast.fitting import trim_grid
from beast.physicsmodel.grid import FileSEDGrid  

if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("trimfile", 
                        help="file with modelgrid, astfiles, obsfiles to use")
    args = parser.parse_args()

    start_time = time.clock()

    f = open(args.trimfile, 'r')
    file_lines = list(f)

    modelfile = file_lines[0].rstrip()
    print('Reading the model grid files = ', modelfile)

    # get the modesedgrid on which to generate the noisemodel  
    modelsedgrid = FileSEDGrid(modelfile)  

    new_time = time.clock()
    print('time to read: ',(new_time - start_time)/60., ' min')

    ext_brick = ''
    
    old_noisefile = ''
    for k in range(1,len(file_lines)):
        line = file_lines[k]
        s1pos = line.find(' ')

        source_density = line[0:s1pos].rstrip()
        sub_source_density = line[s1pos+1:len(line)].rstrip()
    
        basename = '14675_LMC-5665ne-12232.gst'
        obsfile = 'data/' + basename + '_with_sourceden_goodbands' \
                  + '_SD_' + source_density.replace('_','-') + \
                  '_sub' + sub_source_density + '.fits'

        astfile = datamodel.astfile

        noisefile = "%s/%s"%(datamodel.project,datamodel.project) \
                    + '_noisemodel.hd5'

        stats_filebase = "%s/%s"%(datamodel.project,datamodel.project) \
                         + '_sd' + source_density.replace('_','-') \
                         + '_sub' + sub_source_density 
        sed_trimname = stats_filebase + '_sed_trim.grid.hd5'
        noisemodel_trimname = stats_filebase + '_noisemodel_trim.hd5'

        print('working on ' + sed_trimname)
        
        start_time = time.clock()

        if noisefile == old_noisefile:
            print('not reading noisefile - same as last')
            #print(noisefile)
            #print(astfile)
        else:
            print('reading noisefile/astfile')
            # read in the noise model
            noisemodel_vals = noisemodel.get_noisemodelcat(noisefile)
            old_noisefile = noisefile

        # read in the observed data
        print('getting the observed data')
        obsdata = datamodel.get_obscat(obsfile, modelsedgrid.filters)

        # trim the model sedgrid
        #   set n_detected = 0 to disable the trimming of models based on 
        #      the ASTs (e.g. extrapolations are ok)
        #   this is needed as the ASTs in the NIR bands do not go faint enough
        trim_grid.trim_models(modelsedgrid, noisemodel_vals, 
                              obsdata, sed_trimname, noisemodel_trimname, 
                              sigma_fac=3.)

        new_time = time.clock()
        print('time to trim: ',(new_time - start_time)/60., ' min')
