#
# run the orbit code using the namelist version
#
# essential parameters are controlled using a parameter file
#
import run_command as rc
import shutil as SU
import os
import sys

# argumen parsing
import argparse as AG


from LT import parameterfile as PF

# Finding and deleting orb directory. This step is necessary to avoid an error which will occur if run_orbit_nml.py is not fully executed.

if os.path.exists('orb') == True:
	os.remove('orb')
	print "./orb directory has been removed."

parser = AG.ArgumentParser()
parser.add_argument("control_file", nargs = '?', help="Control file ", default = 'control.data')
args = parser.parse_args()

# open control file
p_file = args.control_file

pf = PF.pfile(p_file)

machine = pf.get_value('machine')

input_root = pf.get_value('input_root')

input_dir_ext = pf.get_value('input_dir_ext')

input_dir = input_root + machine + input_dir_ext + '/'

input_file = pf.get_value('input_file')

input_name = input_dir + input_file

input_ext = pf.get_value('nml_ext')

input_file_name = input_name + input_ext

output_root = pf.get_value('output_root')

output_dir_ext = pf.get_value('output_dir_ext')

output_dir = output_root + machine + output_dir_ext + input_file


# create output directory (if necessary)


try:
    os.mkdir(output_dir)
except:
    msg = sys.exc_info()[1]
    if msg.errno == 17 :
        print output_dir, " exists, will use it "
    else:
        print "there is a problem : ", msg
        sys.exit()
# done


print '--------------------------------------------------------------'
print ' Calculation using : ', input_file_name
print 'ALL output in : ', output_dir
print '--------------------------------------------------------------'
# create a symbolic link called ./orb where the default track output of
# orbit goes
try:
    os.symlink(output_dir, "./orb")
except:
    msg = sys.exc_info()[1]
    print "problem with link : ", msg
    sys.exit()

#write the command file


o = open('orbit_input','w')
o.write( input_file_name + '\n')
o.write( 'y\n')
o.close()

# run the code
command = pf.get_value('run_command')
input = open('orbit_input')

ret = rc.run_command(command,'orbit_output', 'orbit_error', stdin = input)

errors = open('orbit_error').readlines()

for l in errors:
    if l.find('Fortran runtime error:') >=0:
        print l
        sys.exit(-1)
    if l.find('No such file or directory') >0:
        print l
        sys.exit(-1)

if ret != 0:
    print 'An error occurred while running orbit, check error file!'
    sys.exit(-1)
    

try:
    # copy the important data file ro the standard  output directory    
    SU.copy('orbits.data', output_dir + '/orbits.data')
    SU.copy('orbit_output', output_dir + '/orbit_output')
    SU.copy('orbit_input', output_dir + '/orbit_input')
    SU.copy('orbit_error', output_dir + '/orbit_error')
    #SU.copy('collimator.data',output_dir + '/collimator.data')
except:
    print '--------------------------------------------------------------'
    print 'problem running orbit'
    print '--------------------------------------------------------------'
    sys.exit(-1)

try:
    print output_dir+'/flux.data'
    SU.copy(input_file_name, output_dir + '/.')
    SU.copy('flux.data', output_dir + '/flux.data')
    SU.copy('flux_limit.data', output_dir + '/flux_limit.data')
    SU.copy('limiter_drawing.data', output_dir + '/limiter_drawing.data')
    SU.copy('fort.8', output_dir + '/fort.8')
except:
    print '--------------------------------------------------------------'
    print 'problem with flux data (orbit)'
    print '--------------------------------------------------------------'
    sys.exit(-1)

#remove temp. files
os.remove('./orbits.data')	
os.remove('./orbit_output')
os.remove('./orbit_input')
os.remove('./orbit_error')
#os.remove('./collimator.data')
os.remove('./orb')

os.remove('./flux.data')
os.remove('./flux_limit.data')
os.remove('./limiter_drawing.data')
os.remove('./fort.8')

#
