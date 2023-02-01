import os
import numpy

from IGMF_analysis.Constants import j1943
# pylint: disable=invalid-name


def _unite_small_files(n_files, file_path, prefix, B_rms, extension='.txt',\
     suffix='_complete', overwrite=False):
    '''Reads all small files with the same characteristics from a folder,
    Joins them into a single larger file, writes it in the same folder and
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
    file_name = (f'{prefix}{B_rms}{extension}')
    large_file = numpy.loadtxt(os.path.join(file_path, file_name))
    length = len(numpy.unique(large_file[11]))
    for i in range(1, int(n_files-1)):
        file_name = (f'{prefix}{i}{B_rms}{extension}')
        if os.path.exists(filename):
            small_file = numpy.loadtxt(os.path.join(file_path, file_name))
            length += len(numpy.unique(small_file[11]))
            if small_file.shape[0] > 0:
                if large_file.ndim == 1:
                    large_file = large_file[:, numpy.newaxis]
                if small_file.ndim == 1:
                    small_file = small_file[:, numpy.newaxis]
                large_file = numpy.concatenate(
                    (large_file, small_file), axis=1)
            else:
                print(f'{file_name} is empty, skipping...')
        else:
            print(f'Warning: {file_name} does not exist')
    print(f'A total of {length} primary photons generated the cascade that\
         hit the detector')
    print(f'Writing output file to {large_file}...')
    return large_file_path


def load_files(n_files, file_path, prefix, B_rms, extension='.txt',
               suffix='_complete'):
    ''' Designed to be chained with the output of unite_small_files, or
    directly with a single large file.
    '''
    return (numpy.loadtxt(_unite_small_files(n_files, file_path, prefix,
            B_rms, extension='.txt', suffix=suffix)))


class simulation:
    '''Class containing all the methods for a CRPropa simulation

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

    def __init__(self, n_files, file_path, prefix, Brms): 
        ''' This has to improve. Create a dict and loop over its entries
        '''
        table = load_files(n_files, file_path, prefix, Brms)
        self.photon_mask = numpy.abs(table[3]) == 22
        self.cascade_mask = ~(table[2] == table[11])
        self.cascade_photon = self.photon_mask*self.cascade_mask
        self.px = table[8]
        self.py = table[9]
        self.pz = table[10]
        self.px0 = table[17]
        self.py0 = table[18]
        self.pz0 = table[19]
        self.energy = table[4]*1e9
        self.idx = numpy.abs(table[3]) == 22
        self.x0 = table[14]
        self.y0 = table[15]
        self.z0 = table[16]
        self.x1 = table[23]
        self.y1 = table[24]
        self.z1 = table[25]
        self.x = table[5]
        self.y = table[6]
        self.z = table[7]

    def dermer_cut(self, angle):
        ''' Utility function to cut the photons based on the angle as
        defined in Dermer et al (2011) paper
        https://ui.adsabs.harvard.edu/abs/2011ApJ...733L..21D/abstract.
        '''
        # Create the vectors to get the angle delta
        delta_0 = numpy.transpose((self.px, self.py, self.pz))
        delta_1 = numpy.transpose((self.px0, self.py0, self.pz0))
        # delta from the dot product:
        # clip to limit the extrema in [-1,1], not mathematically required
        # but floating point approximations screw up the process
        # multiply to make a parallel calculation element-wise
        # sum gieves the element-wise summation on the given axis
        delta = numpy.array(numpy.arccos(numpy.clip(numpy.multiply\
                            (delta_0, delta_1).sum(1)), -1, 1))
        lambda_xx = numpy.sqrt((self.x1-self.x0)**2 + (self.z1-self.z0)**2
                               + (self.y1-self.y0)**2)
        # The following is stored in the class for future usage (e.g.: plotting)
        self.dermer_theta = numpy.arcsin(
            (lambda_xx/j1943.D) * numpy.cos(numpy.pi/2 - delta))
        # Define the hit condition to filter photons based on the angle
        hit = numpy.sqrt(self.dermer_theta*(180./numpy.pi)**2) < angle
        return (hit)

    def momentum_cut(self, angle=0.3):
        ''' Utility function to cut the photons based on the angle as
        derived from the momentum of the incoming photon.

        LAT angle to be defined in (fractional) degrees and is 0.3 by default
        '''
        p_norm = self.px**2 + self.py**2 + self.pz**2
        x_arr = (self.px/p_norm)
        y_arr = (self.py/p_norm)
        # This other one has been entirely rewritten: the angle is the
        # weight of the projection w.r.t. the entire vector amplitude,
        # which is 1 because the vector is normalized. We don't distinguish
        # between x and y for the angular opening, so we use them both.
        # angle should be 0.3 for the LAT
        hit = numpy.sqrt(numpy.arcsin(x_arr)**2 + numpy.arcsin(y_arr)**2)\
            * 180/numpy.pi < angle
        return (hit)
