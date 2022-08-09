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