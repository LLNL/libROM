import numpy
import h5py
import argparse

parser = argparse.ArgumentParser(
            prog='Basis data format converter',
            description='This python script converts the out-dated dataset names of '
                        'a basis hdf5 file. 6 digits of integer following '
                        'the name of all datasets will be removed.')
parser.add_argument('filename', type=str,
                    help='the basis hdf5 file name to update the dataset names.')           # positional argument

if __name__ == "__main__":
    args = parser.parse_args()
    print('Filename: %s' % args.filename)
    with h5py.File(args.filename, 'r+') as f:
        for key in f.keys():
            if ((key[-7] == '_') and (key[-6:].isdigit())):
                print('Detected an outdated dataset name: %s' % key)
                if (int(key[-6:]) != 0):
                    print('libROM only supports single time interval, skipping this dataset.')
                    continue
                newname = key[:-7]
                print('Converting to: %s' % newname)
                f[newname] = f[key]
                del f[key]


