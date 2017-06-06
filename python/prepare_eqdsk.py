# Prepare eqdsk file for orbit
#
# shell utilities
#
# run this script in the orbit3_gfortran directory
#
import shutil as SH
import run_command as rc
import glob as G
import os.path as P

# fortran program to be used
edit_eqdsk = './fortran/edit_eqdsk > eqdsk_new.dat'

# source directory
src_dir = './MAST_efit/efit_jstorrs/26789/'
# use all files that start with g
src_pattern = 'g*'

# get all file names
eqdsk_files = G.glob(src_dir + src_pattern)


dest_dir = './MAST_efit/efit_jstorrs/'

# now loop over all file_names
for ef in eqdsk_files:
    # copy the file to eqdsk_old.dat
    SH.copy(ef, 'eqdsk_old.dat')
    # run the fortran editor
    rc.run_command( edit_eqdsk, 'edit_eqdsk.err', 'edit_eqdsk.out')
    # copy the result to its final location
    fname = P.split(ef)[-1]
    SH.copy('eqdsk_new.dat', dest_dir + fname + '.dat')
#
# all done


    