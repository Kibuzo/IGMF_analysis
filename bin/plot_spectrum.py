import logging
from argparse import ArgumentParser
import os

import numpy
import matplotlib.pyplot as plt

from IO import spectrum
from math_.base import cutoff
from structure import IGMF_DATA

# pylint: disable=logging-fstring-interpolation
# pylint: disable=invalid-name

__description__ = \
"""Utility to plot the spectra from our unified standard fits file.
All files should have the same format, independently of the source that
produced them

As of now, this is just a structure for the file.
"""
DATA_SET = ['fermi', 'veritas', 'hess', 'fit', 'cascade', 'all']
PARSER = ArgumentParser(description=__description__)
PARSER.add_argument('--include', nargs = '*', choices=DATA_SET, default='all', \
    help='Plots to include and overlay')

def plot_cutoff (energy = [7,4060], N0=3.83e-14, index=1.5 , ecut=2080):
    ''' Pass energy grid, which will then be resampled to create a smoother
    source function plot.
    The function is the expcutoff defined in math.base, and the default
    parameters are the ones defined in the paper.
    '''
    e_grid = numpy.geomspace(min(energy), max(energy), 1000)
    f_grid = cutoff(e_grid, N0, index, ecut)
    plt.plot(e_grid, f_grid, ls = '--', color='tab:orange', \
        label = 'Exponential cutoff fit')


def make_plot(plot_list):
    ''' Main plotting app. It should plot from fits files since we have defined
    universal format for that reason. Specific data sets point to specific files
    path so that by choosing a name we are actually implying a path.

    Warning, this is not designed for flexibility, as it only provides a quick
    interface to plot static data, including the fits. If you want to test stuff
    do it manually importing the libraries instead.
    '''
    spec_name = {'fermi': 'Fermi.fits', 'hess': 'Hess.fits', \
         'veritas': 'Veritas.fits','cascade': 'cascade.fits', \
         'fit': 'fit.fits'}
    spec_symbols = {'fermi': '.', 'hess': 'v', 'veritas': 'x'}
    if plot_list == 'all':
        plot_list = DATA_SET[:-1]
    for plot in plot_list:
        logging.info(f'Attempting to plot {plot}...')
        if plot == 'fit':
            plot_cutoff()
        if plot == 'cascade':
            pass
        else:
            try:
                spec = spectrum.from_fits(os.path.join(IGMF_DATA, \
                    spec_name[plot]))
                spec.plot(marker = spec_symbols[plot], color = 'tab:blue', \
                    label = f'de-absorbed data points ({plot.upper()})')
            except:
                logging.warning(f'Plotting of {plot} failed, maybe the fits '\
                    'file is not in the data folder?')
    plt.title ('HESS J1943+213 intrinsic spectrum')
    plt.legend()
    plt.show()

def main():
    """main() entry point.
    """
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    args=PARSER.parse_args().__dict__
    plt.figure('Spectrum')
    make_plot (args['include'])

if __name__ == '__main__':
    main()
