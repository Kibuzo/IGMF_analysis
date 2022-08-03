import numpy
import os

def unite_small_files(n_files, file_path, prefix, B_rms, suffix = '.txt'):
    '''Reads all small files with the same characteristics from a folder,
    Joins them into a single larger file and returns the corresponding
    array.

    This sucks badly and over time recieved several patches decreasing
    the readability even more. Now it's unreadable, but works.
    Consider destroying it entirely and rewriting it.
    '''
    large_file=[] #This needs to exist in case first file is empty
    file_name = (f'{pefix}{B_rms}{suffix}')
    large_file=numpy.loadtxt(os.path.join(file_path, file_name))
    length=len(numpy.unique(a[11]))
    for i in range (1,int(n_files-1)):
        file_name=(f'{prefix}{i}{B_rms}{suffix}')
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
    return (large_file)

def load_large_file(file_path, prefix, Brms, suffix='.txt'):
    file_name = (f'{prefix}{Brms}{suffix}')
    return (numpy.loadtxt(os.path.join(file_path, file_name)))