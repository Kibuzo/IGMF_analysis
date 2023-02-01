import logging
import IO.spectrum
from argparse import ArgumentParser
import matplotlib.pyplot as plt
__description__ = \
"""Utility to plot the spectra from our unified standard fits file.
All files should have the same format, independently of the source that
produced them

As of now, this is just a structure for the file.
"""
DATA_SET = ['fermi', 'veritas', 'hess', 'fit', 'cascade','all']
PARSER = ArgumentParser(description=__description__)
PARSER.add_argument('--include', nargs = '*', choices=DATA_SET, default='all', \
    help='Plots to include and overlay')



def plot(plot_list):
    ''' Main plotting app. It should plot from fits files since we have defined
    universal format for that reason. Specific data sets point to specific files
    so that by choosing a name we are actually implying a path.
    '''
    for plot in plot_list:
        logging.info(f'plotting {plot}...')
    

def main():
    """main() entry point.
    """
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    args=PARSER.parse_args().__dict__
    plt.figure('Spectrum')
    plot (args['include'])

if __name__ == '__main__':
    main()