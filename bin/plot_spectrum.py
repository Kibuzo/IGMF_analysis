import logging
import os
from IO import spectrum
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from structure import IGMF_DATA


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



def make_plot(plot_list):
    ''' Main plotting app. It should plot from fits files since we have defined
    universal format for that reason. Specific data sets point to specific files
    so that by choosing a name we are actually implying a path.
    '''
    spec_name = {'fermi': 'Fermi.fits', 'hess': 'Hess.fits', 'veritas': 'Veritas.fits', \
        'cascade': 'cascade.fits', 'fit': 'fit.fits'}
    if plot_list == 'all':
        plot_list = DATA_SET[:-1]
    for plot in plot_list:
        logging.info(f'Attempting to plot {plot}...')
        try:
            spec = spectrum.from_fits(os.path.join(IGMF_DATA, spec_name[plot]))
            spec.plot(marker = '.', color = 'tab:blue')
        except:
            logging.warning(f'Plotting of {plot} failed, maybe the fits file is\
                not in the data folder?')
            pass
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