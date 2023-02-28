import logging
import os
import numpy
import matplotlib.pyplot as plt

from Constants.J1943 import D
from math_.base import xyoffset_to_r, radians_to_degree, degree_to_radians
# pylint: disable=invalid-name

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

CRPROPA_DICT_NAMES = ['D', 'z',	'SN', 'ID',	'E', 'X', 'Y', 'Z',	'Px', 'Py',	'Pz',\
    'SN0', 'ID0', 'E0', 'X0', 'Y0', 'Z0', 'P0x', 'P0y', 'P0z', 'SN1', 'ID1', 'E1',\
    'X1', 'Y1', 'Z1', 'P1x', 'P1y', 'P1z', 'W']
CRPROPA_NUM = len(CRPROPA_DICT_NAMES)

def _unite_small_files(n_files, file_path, prefix, B_rms, extension='.txt',\
     suffix='_complete', overwrite=False):
    '''Reads all small files with the same characteristics from a folder,
    joins them into a single larger file, writes it in the same folder and
    returns the corresponding file path and name.

    This sucks badly and over time recieved several patches decreasing
    the readability even more. Now it's unreadable, but works.
    Consider destroying it entirely and rewriting it.
    Keep in mind that its design was thought as necessary. The chain of
    'ifs'  is due to the fact that empty files might occour and disrupt
    the reading pipeline with arrays that have wrong dimension, so an
    accurate refactoring should account for this possibility.

    The way this was originally designed to work was this: 
        - Start by simulating a small subset of photons
        - Excise the useless columns and place everything in an ascii file
        with no labels at all (numpy's savetxt was the powerhorse here...)
        - copy all the files from the cluster to the local folder
        - Recover the original formatting and unite everything in a single 
        large file still with no formatting
    I should try at least to document the structure and provide public, 
    cleanedand possibly sort of human-readable files for all simulations 
    that i have performed to get rid of local copies, and this function 
    could be the place to start..
    A temporary file archive on which to run this script can be found at
    https://drive.google.com/file/d/1qPSTzlKU7pKz00GVkVVYc1qIPqSN1b7G/view?usp=share_link
    '''
    large_file_name = (f'{prefix}{B_rms}{suffix}{extension}')
    large_file_path = os.path.join(file_path, large_file_name)
    if not overwrite and os.path.exists(large_file_path):
        return large_file_path
    large_file = []  # This needs to exist in case first file is empty
    file_name = (f'{prefix}{0}{B_rms}{extension}')
    large_file = numpy.loadtxt(os.path.join(file_path, file_name))
    length = len(numpy.unique(large_file[11]))
    for i in range(1, int(n_files-1)):
        #skip reading the first one cause it's read outside the loop
        file_name = (f'{prefix}{i}{B_rms}{extension}')
        full_path = os.path.join(file_path, file_name)
        if os.path.exists(full_path):
            logging.info(f'reading from {full_path}...')
            small_file = numpy.loadtxt(os.path.join(file_path, file_name))
            # Horrible hack to skip empty files
            if os.stat(full_path).st_size == 30:
                logging.warning(f'file {small_file} appears to be empty,'
                                'skipping...')
                continue
            length += len(numpy.unique(small_file[11]))
            if small_file.shape[0] > 0:
                if large_file.ndim == 1:
                    large_file = large_file[:, numpy.newaxis]
                if small_file.ndim == 1:
                    small_file = small_file[:, numpy.newaxis]
                large_file = numpy.concatenate(
                    (large_file, small_file), axis=1)
            else:
                logging.info(f'{file_name} is empty, skipping...')
        else:
            logging.warning(f'{file_name} does not exist')
    logging.info(f'A total of {length} primary photons generated the '\
        'cascade that hit the detector')
    logging.info(f'Writing output file to {large_file_path}...')
    numpy.savetxt(large_file_path, large_file)
    return (large_file_path)


def load_files(n_files, file_path, prefix, B_rms, extension='.txt',
               suffix='_complete', overwrite = False):
    ''' Designed to be chained with the output of unite_small_files, or
    directly with a single large file.
    '''
    return (numpy.loadtxt(_unite_small_files(n_files, file_path, prefix,
            B_rms, extension=extension, suffix=suffix, overwrite = overwrite)))


class simulation:
    '''Class containing all the methods for a CRPropa simulation.
    The CRPropa default 3D simulation column output is as follows
    
    D	z	SN	ID	E	X	Y	Z	Px	Py	Pz	SN0	ID0	E0	X0	Y0	Z0
    P0x	P0y	P0z	SN1	ID1	E1	X1	Y1	Z1	P1x	P1y	P1z	W

    D             Trajectory length [1 Mpc]
    z             Redshift
    SN/SN0/SN1    Serial number. Unique (within this run) id of the particle.
    ID/ID0/ID1    Particle type (PDG MC numbering scheme)
    E/E0/E1       Energy [1 EeV]
    X/X0/X1...    Position [1 Mpc]
    Px/P0x/P1x... Heading (unit vector of momentum)
    W             Weights
    no index = current, 0 = at source, 1 = at point of creation
    '''

    def __init__(self, n_files, file_path, prefix, Brms, overwrite = False): 
        ''' This has to improve. Create a dict and loop over its entries
        '''
        table = load_files(n_files, file_path, prefix, Brms, overwrite=overwrite)
        self.data = {CRPROPA_DICT_NAMES[i]: table[i] for i in range (CRPROPA_NUM)}
        #From EeV to GeV
        self.data['E'] = self.data['E']*1e9
        self.data['E1'] = self.data['E1']*1e9
        self.data['E0'] = self.data['E0']*1e9
        logging.warning('energies are now in GeV')
        self.photon_mask = numpy.abs(table[3]) == 22
        self.cascade_mask = ~(table[2] == table[11])
        self.cascade_photon = self.photon_mask*self.cascade_mask
        # Demer 11 variables
        # Create the vectors to get the angle delta
        delta_0 = numpy.transpose((self.data['Px'], self.data['Py'], self.data['Pz']))
        delta_1 = numpy.transpose((self.data['P0x'], self.data['P0y'], self.data['P0z']))
        # delta from the dot product:
        # clip to limit the extrema in [-1,1], not mathematically required
        # but floating point approximations screw up the process
        # multiply to make a parallel calculation element-wise
        # sum gives the element-wise summation on the given axis
        delta = numpy.array(numpy.arccos(numpy.clip(numpy.multiply 
                            (delta_0, delta_1).sum(1), -1, 1)))
        lambda_xx = numpy.sqrt((self.data['X1']-self.data['X0'])**2 + (self.data['Z1']-self.data['Z0'])**2
                               + (self.data['Y1']-self.data['Y0'])**2)
        self.dermer_theta = numpy.arcsin(
            (lambda_xx/D) * numpy.cos(numpy.pi/2 - delta))
        nphot = len(numpy.unique(self.data['SN0']))
        logging.info(f'The simulations contains a total of {nphot} source photons')


    
    def photon_mask(self):
        '''
        '''
        pass

    def write_fits(self, file_path):
        ''' Interface for saving the simulation in the form of a fits,
        compatible with those of the telescope's spectral analysis.
        The fit should have a spectrum column and possibly some other
        observable  that we would like to plot after the selection cut 
        performed by our processing
        '''
        pass


    def dermer_cut(self, angle):
        ''' Utility function to cut the photons based on the angle as
        defined in Dermer et al. (2011) paper
        https://ui.adsabs.harvard.edu/abs/2011ApJ...733L..21D/abstract.
        '''
        # Define the hit condition to filter photons based on the angle
        self.dermer_cut = numpy.sqrt(self.dermer_theta*(180./numpy.pi)**2)\
                < angle

    def momentum_cut(self, angle=0.3):
        ''' Utility function to cut the photons based on the angle as
        derived from the momentum of the incoming photon.

        LAT angle to be defined in (fractional) degrees and is 0.3 by default

        Note that this is not meant to be a selection of "Good" photons in the
        case in which we are plotting a map
        '''
        p_norm = self.data['Px']**2 + self.data['Py']**2 + self.data['Pz']**2
        x_arr = (self.data['Px']/p_norm)
        y_arr = (self.data['Py']/p_norm)
        # This other one has been entirely rewritten: the angle is the
        # weight of the projection w.r.t. the entire vector amplitude,
        # which is 1 because the vector is normalized. We don't distinguish
        # between x and y for the angular opening, so we use them both.
        # angle should be 0.3 for the LAT

        #hit = numpy.sqrt(numpy.arcsin(x_arr)**2 + numpy.arcsin(y_arr)**2)\
        #    * 180/numpy.pi < angle

        self.momentum_cut = radians_to_degree(xyoffset_to_r(numpy.arctan(x_arr), \
            numpy.arctan(y_arr)**2)) < angle
    
    def coordinate_cut(self, half_angle = 2.5):
        ''' Select only photons hitting the (virtual) detector. This is where
        you want to start from if you want to create a map
        '''
        self.hit = radians_to_degree(xyoffset_to_r(self.data['X'], self.data['Y']))\
                 < half_angle
    
    def plot_radial_profile (self, max_rad = 2, dermer = False, **kwargs):
        ''' Create the radial profile of the incoming photons, defined a-la-dermer
        or not with the appropriate flag. max_rad is needed to have consistent
        binning between the two versions.
        '''
        if dermer:
            plt.figure('Radial profile (Dermer plot)')
            radial = (radians_to_degree(self.dermer_theta))
            plt.hist(radial[radial<max_rad], weights=1./radial[radial<max_rad],
                     **kwargs)
            plt.yscale('log')
            plt.ylabel('Area-weighted counts')
            plt.xlabel('Angular distance [degrees]')
        else:
            theta_x = radians_to_degree(self.data['Px'])
            theta_y = radians_to_degree(self.data['Py'])
            theta_r = numpy.sqrt(theta_x**2+theta_y**2)
            plt.figure ('Radial profile')
            plt.hist(theta_r[theta_r<max_rad], weights=1./theta_r[theta_r<max_rad],
                     **kwargs)
            plt.yscale('log')
            plt.ylabel('Area-weighted counts')
            plt.xlabel('Angular distance taken from beam center[degrees]')

    def plot_map (self, selection = 'momentum'):
        ''' 2d theta² plot of the map of incoming photons.
        '''
        theta_x = radians_to_degree(self.data['Px'])
        theta_y = radians_to_degree(self.data['Py'])
        try:
            plt.figure('momentum cut')
            plt.hist2d (theta_x[self.hit], theta_y[self.hit], bins=360)
            plt.xlabel(r'$\theta_x$')
            plt.ylabel(r'$\theta_y$')
            plt.title('Arrival direction of photons (momentum cut)')
        except:
            logging.warning ('No momentum cut performed, run'\
                            'simulation.momentum_cut() before running this')
        try:
            plt.figure('dermer cut')
            plt.hist2d (theta_x[self.hit], theta_y[self.hit], bins=360)
            plt.xlabel(r'$\theta_x$')
            plt.ylabel(r'$\theta_y$')
            plt.title('Arrival direction of photons (D11 cut)')
        except:
            logging.warning ('No Dermer cut performed, run'\
                            'simulation.momentum_cut() before running this')
        

    def plot_raytrace(self):
        ''' 2d plot of a slice of the space that displays the interaction
        point and the direction of electromagnetic particles.
        '''
        pass
