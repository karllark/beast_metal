#!/usr/bin/env python
""" 
Script to plot the variation in the noise model versus source density
"""
import argparse

import numpy as np
import tables
import matplotlib.pyplot as plt

import datamodel
import beast.observationmodel.noisemodel.generic_noisemodel as noisemodel 
from beast.observationmodel.vega import Vega
from beast.plotting.beastplotlib import initialize_parser

if __name__ == '__main__':

    # commandline parser
    parser = initialize_parser()
    parser.add_argument("astfile", type=str,
                        help="Artifical Star Test (AST) FITS file")
    parser.add_argument("--mag_out", action="store_true",
                        help="Use output mags instead of rate [flux]")
    parser.add_argument("--model_noise_file", type=str, default=None,
                        help="Model noise model file")
    args = parser.parse_args()
    astfile = args.astfile
    basename = astfile.replace('.fits','')

    filters = datamodel.basefilters
    filters_tag = datamodel.basefilters
    n_filters = len(datamodel.filters)

    show_indiv = True

    with Vega() as v:
        _, vega_flux, _ = v.getFlux(datamodel.filters)

    fig, ax = plt.subplots(n_filters, sharex=True, figsize=(10,13))
    
    blabels = ['all']

    colsym = ['co','go','ro','bo','yo']
    col = ['c','g','r','b','y']

    for k in range(len(blabels)):
        blabel = blabels[k]
        model = noisemodel.Generic_ToothPick_Noisemodel(astfile, filters)

        #model.fit_bins(nbins=30, completeness_mag_cut=-10)
        
        for i in range(n_filters):
            mag_in = model.data[filters[i] + '_in']

            if args.mag_out:
                completeness_mag_cut = 70.
            else:
                # change mappings for out column to the rate column
                model.data.set_alias(filters[i] + '_out',
                                     filters[i] + '_RATE')
                completeness_mag_cut = -10.
            mag_out = model.data[filters[i] + '_out']

            if show_indiv:
                flux_in_all = 10 ** (-0.4*mag_in)
                if args.mag_out:
                    flux_out_all = 10 ** (-0.4*mag_out)
                    bad_indxs,= np.where(mag_out >= completeness_mag_cut)
                    flux_out_all[bad_indxs] = 0.0
                else:
                    flux_out_all = mag_out

                indxs, = np.where(flux_out_all != 0.0)
                ax[i].plot(vega_flux[i]*flux_in_all[indxs],
                           flux_out_all[indxs]/flux_in_all[indxs] - 1.0,'ko',
                           markerfacecolor='none',markeredgecolor='k')

                indxs, = np.where(abs(flux_out_all) == 0.0)
                ax[i].plot(vega_flux[i]*flux_in_all[indxs],
                           flux_out_all[indxs]/flux_in_all[indxs] - 1.0,
                           'ko')
            
            d = model._compute_sigma_bins(mag_in, mag_out, nbins=50, 
                                    completeness_mag_cut=completeness_mag_cut)

            if i == 0:
                ax[i].plot(vega_flux[i]*d['FLUX_IN'],d['FLUX_BIAS']
                           /d['FLUX_IN'],colsym[k], label=blabel + ' $\mu$')
            else:
                ax[i].plot(vega_flux[i]*d['FLUX_IN'],
                           d['FLUX_BIAS']/d['FLUX_IN'],
                           colsym[k])

            ax[i].plot(vega_flux[i]*d['FLUX_IN'],d['FLUX_BIAS']
                       /d['FLUX_IN'],colsym[k].replace('o','-'))

            if i == 1:
                ax[i].plot(vega_flux[i]*d['FLUX_IN'],(d['FLUX_STD']
                                                      /d['FLUX_IN']),colsym[k],
                           markerfacecolor='none',markeredgecolor=col[k], 
                           label=blabel + ' $\sigma$')
            else:
                ax[i].plot(vega_flux[i]*d['FLUX_IN'],(d['FLUX_STD']
                                                      /d['FLUX_IN']),colsym[k],
                           markerfacecolor='none',markeredgecolor=col[k])

            ax[i].plot(vega_flux[i]*d['FLUX_IN'],(d['FLUX_STD']
                                                  /d['FLUX_IN']),
                       colsym[k].replace('o','-'),
                       markerfacecolor='none',markeredgecolor=col[k])

            ax[i].text(1e-16,2.5,filters_tag[i])

    for i in range(n_filters):
        ax[i].set_yscale('linear')
        #ax[i].set_ylabel('($\mu$ or $\sigma$)/F(band)') 
        ax[i].set_ylabel('$\mu$/F(band)') 
        #ax[i].set_ylabel('Sigma/Flux(Input)')
        #ax[i].set_ylim(-0.5,2.0)
        ax[i].set_ylim(-3,3)
        ax[i].set_xscale('log')
        ax[i].set_xlim(1e-25,2e-14)

    ax[n_filters-1].set_xlabel('F(band) [ergs cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]')
    ax[0].legend()
    ax[1].legend()

    plt.tight_layout()
    fig.subplots_adjust(hspace=0.05)

    outbase = basename + '_ast_model'
    if args.mag_out:
        outbase += '_mag_out'

    if args.savefig:
        fig.savefig('{}.{}'.format(outbase, args.savefig))
    else:
        plt.show()

