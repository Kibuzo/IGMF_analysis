import numpy
import os

def _unite_small_files(n_files, file_path, prefix, B_rms, extension =\
    '.txt',suffix='_complete', overwrite=False):
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
    '''
    large_file_name = (f'{prefix}{B_rms}{suffix}{extension}')
    large_file_path = os.path.join(file_path, large_file_name)
    if not overwrite and os.path.exists(large_file):
        return large_file_path
    large_file=[] #This needs to exist in case first file is empty
    file_name = (f'{pefix}{B_rms}{extension}')
    large_file=numpy.loadtxt(os.path.join(file_path, file_name))
    length=len(numpy.unique(a[11]))
    for i in range (1,int(n_files-1)):
        file_name=(f'{prefix}{i}{B_rms}{extension}')
        if os.path.exists(filename):
            small_file=numpy.loadtxt(os.path.join(file_path, file_name))
            length+=len(numpy.unique(s[11]))
            if small_file.shape[0]>0:
                if (large_file.ndim==1):
                    large_file=large_file[:, numpy.newaxis]
                if (small_file.ndim==1):
                    small_file=small_file[:, numpy.newaxis]
                large_file=numpy.concatenate((large_file,small_file),axis=1)
            else:
                print (f'{file_name} is empty, skipping...')
        else:
            print (f'Warning: {file_name} does not exist')
    print(f'A total of {length} primary photons generated the cascade that\
         hit the detector')
    print(f'Writing output file to {large_file}...') 
    return (large_file_path)

def load_files(n_files, file_path, prefix, B_rms, extension='.txt',
               suffix='_complete'):
    ''' Designed to be chained with the output of unite_small_files, or
    directly with a single large file.
    '''
    return (numpy.loadtxt(_unite_small_files(n_files, file_path, prefix, B_rms, 
            extension = '.txt',suffix = suffix)))

class simulation:
    '''Class containing all the methods for a CRPropa simulation
    '''
    def __init__(self, file_path, n_files):
        table = load_files(n_files, file_path)
        self.photon_mask = numpy.abs(table[3]) == 22
        self.cascade_mask = ~(table[2]==table[11])
        self.cascade_photon = self.photon_mask*self.cascade_mask
        self.px = table[8]
        self.py = table[9]
        self.pz = table[10]
        self.px0 = table[17]
        self.py0 = table[18]
        self.pz0 = table[19]
        self.energy = table[4]*1e9
        self.idx = numpy.abs(table[3]) == 22
        pnorm = numpy.sqrt(self.px**2+self.py**2+self.pz**2)
        self.x0 = table[14]
        self.yo = table[15]
        self.z0 = table[16]
        self.x1 = table[23]
        self.y1 = table[24]
        self.z1 = table[25]
        self.x = table[5]
        self.y = table[6]
        self.z = table[7]



    def neronov_cut(self, angle):
        ''' Utility function to cut the photons based on the angle as
        defined in Neronov paper.
        '''
        delta = numpy.transpose((self.px, self.py, self.pz))
        delta_initial = numpy.transpose((self.px0, self,py0, self.pz0))
        #Vectors are already normalized in CRPropa
        dot_product = numpy.array(numpy.arccos(numpy.clip(numpy.multiply(delta,
                      delta_initial).sum(1)), -1, 1))
        phi = self.y1 / numpy.sqrt(self.y1**2+self.x1**2)
        cos_phi = numpy.sign(self.x1)*numpy.cos(numpy.array(phi))
        sin_phi = numpy.sign(y1)*numpy.sqrt(1-cos_phi**2)
        lambda_xx = numpy.sqrt((self.x1-self.x0)**2 + (self.z1-self.z0)**2 + (self.y1-self.y0)**2)
        theta = numpy.arcsin((lambda_xx/d) * numpy.cos(numpy.pi/2 - dot_product))
        theta_x = theta*cos_phi
        theta_y = theta*sin_phi
        #Define the hit condition to filter photons based on the angle
        self.hit = numpy.sqrt(theta*(180./numpy.pi)**2)<angle

    def coordinate_cut(self, angle):
        ''' Utility function to cut the photons based on the angle as
        derived from the momentum of the incoming photon.
        '''
        pass



