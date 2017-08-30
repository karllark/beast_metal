#!/usr/bin/env python

"""
Code to setup the xsede slurm inputs for results that are incomplete
  incomplete means 1D pPDF results done with the old linear mass spacing
    or missing full 1D pPDF or lnp files

.. history::
    Written 4jul15 by KDG.

"""

from __future__ import print_function
import os
import glob
import argparse
import tables
#####
import numpy as np
from astropy.table import Table
from astropy.io import fits
#####

if __name__ == '__main__':

    basename = '14675_LMC-5665ne-12232.gst_with_sourceden_goodbands'
    datapath = 'data'
    basepath = 'LMC-5665ne-12232_beast'
    cat_files = np.array(glob.glob(datapath+'/'+basename + '*_sub*.fits'))

    n_cat_files = len(cat_files)
    n_pernode_files = 1

    # setup the subdirectory for the xsede slurm and log files
    job_path = basepath+'/refit_batch_jobs/'
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

    log_path = job_path+'logs/'
    if not os.path.isdir(log_path):
        os.mkdir(log_path)

    pf_open = False
    cur_f = 0
    cur_total_size = 0.0
    j = -1

    #cat_files = cat_files[0:2]

    for i, cat_file in enumerate(cat_files):
        # get the sd number
        bpos = cat_file.find('obscat/')
        dpos = cat_file.find('SD_')
        spos = cat_file.find('sub')
        ppos = cat_file.rfind('.')
        sd_num = cat_file[dpos+3:spos-1]
        sub_num = cat_file[spos+3:ppos]

        # read the stats file and see if this subregion is done yet
        results_path = basepath + 'results/'
        stats_file = results_path+'http_sd'+sd_num+'_sub'+sub_num+'_stats.fits'
        pdf1d_file = results_path+'http_sd'+sd_num+'_sub'+sub_num+'_pdf1d.fits'
        lnp_file = results_path+'http_sd'+sd_num+'_sub'+sub_num+'_lnp.hd5'

        reg_run = False
        run_done = False
        if not os.path.isfile(stats_file):
            reg_run = True
            print('no stats file')
        if not os.path.isfile(pdf1d_file):
            reg_run = True
            print('no pdf1d file')
        if not os.path.isfile(lnp_file):
            reg_run = True
            print('no lnp file')

        # first check if the pdf1d mass spacing is correct
        if not reg_run:
            hdulist = fits.open(pdf1d_file)
            delta1 = hdulist['M_ini'].data[-1,1] - hdulist['M_ini'].data[-1,0]
            if delta1 > 1.0:  # old linear spacing
                print('pdf1d lin mass spacing - full refitting needed')
                old_mass_spacing = True
            else:
                old_mass_spacing = False
                print('pdf1d log mass spacing - ok')

            if old_mass_spacing:
                run_done = False
                reg_run = True

        # now check if the number of results is the same as 
        #    the number of observations
        if not reg_run:
            # get the observed catalog
            obs = Table.read(cat_file)

            # get the fit results catalog
            t = Table.read(stats_file)
            # get the number of stars that have been fit
            indxs, = np.where(t['Pmax'] != 0.0)

            # get the number of entries in the lnp file
            f = tables.openFile(lnp_file, 'r')
            nlnp = f.root._v_nchildren - 2
            f.close()

            print('# obs, stats, lnp = ',len(obs), len(indxs), nlnp)
            if (len(indxs) == len(obs)) & (nlnp == len(obs)):

                # final check, is the pdf1d file correctly populated
                tot_prob = np.sum(hdulist['M_ini'].data, axis=1)
                tindxs, = np.where(tot_prob > 0.0)
                print('# good pdf1d = ', len(tindxs) - 1)
                if len(tindxs) == (len(obs) + 1):
                    run_done = True

        if run_done:
            print(stats_file + ' done')
        else:

            j += 1
            if j%n_pernode_files == 0:
                cur_f += 1

                # close previous files
                if j != 0:
                    pf.close()
                    print('total sed_trim size [Gb] = ', 
                          cur_total_size/(1024.*1024.*1024.))
                    cur_total_size = 0.0

                # open the slurm and param files
                pf_open = True
                joblist_file = job_path+'beast_batch_refit_'+str(cur_f) \
                               +'.joblist'
                pf = open(joblist_file,'w')
                

            ext_str = ''
            if reg_run:
                print(stats_file 
                      + ' does not exist - adding job as a regular fit job (not resume job)')
            else:
                print(stats_file 
                      + ' not done - adding to continue fitting list (' + \
                      str(len(indxs)) + '/' + str(len(t['Pmax'])) + ')')
                ext_str = '-r'

            job_command = './run_beast_production.py -f ' + ext_str + ' ' + \
                          sd_num + ' '+sub_num+' > ' \
                          + log_path+'beast_fit_resume_http' + \
                          '_sd'+sd_num+'_sub'+sub_num+'.log'

            pf.write(job_command+'\n')

    if pf_open:
        pf.close()
