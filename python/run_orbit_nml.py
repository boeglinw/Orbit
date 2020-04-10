#
# run the orbit code using the namelist version
#
import run_command as rc
import shutil as SU
import os
import sys

# Finding and deleting orb directory. This step is necessary to avoid an error which will occur if run_orbit_nml.py is not fully executed.

if os.path.exists('orb') == True:
	os.remove('orb')
	print "./orb directory has been removed."

#machine = 'NSTX'
machine = 'MAST'

input_dir = './' + machine + '_input/'


#input_file = 'MAST_p2'
#input_file = 'ppro_6'
#input_file= 'p_test'

# parameter files for renewal proposal
#input_file = 'ppro_7'
#input_file = 'ppro_7_bundle'
#input_file = 'ppro_7_high_B'
#input_file = 'ppro_7_bundle_high_B'

#input_file = 'p8_upper_position_high_B'
#input_file = 'p8_lower_position_high_B'

#input_file = 'p7_upper_position_high_B'
#input_file = 'p7_mid_position_high_B'
#input_file = 'p7_mid_position'

#input_file = 'p7_lower_position_high_B'
#input_file = 'p7_bundle_lower_position_high_B'

#input_file = 'p6_lower_position'
#input_file = 'ppro_5'
#input_file = 'p_lower_position'
#input_file = 'ppro_8'
#input_file = 'ppro_6'

# MAST calculations
#input_file = 'MAST_p2_c'
input_file = 'MAST_p6'

input_name = input_dir + input_file

input_ext = '.nml'

output_dir = './'+ machine + '_output/nml_orb_'+input_file
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
o.write( input_name + input_ext + '\n')
o.write( 'y\n')
o.close()

# run the code
command = './bin/orbit205'
input = open('orbit_input')
ret = rc.run_command(command,'orbit_output', 'orbit_error', stdin = input)

# copy the important data file ro the standard  output directory
SU.copy('orbits.data', output_dir + '/orbits.data')
SU.copy('orbit_output', output_dir + '/orbit_output')
SU.copy('orbit_input', output_dir + '/orbit_input')
SU.copy('orbit_error', output_dir + '/orbit_error')
#SU.copy('collimator.data',output_dir + '/collimator.data')

os.remove('./orbits.data')	
os.remove('./orbit_output')
os.remove('./orbit_input')
os.remove('./orbit_error')
#os.remove('./collimator.data')
os.remove('./orb')

print output_dir+'/flux.data'
SU.copy('flux.data', output_dir + '/flux.data')
SU.copy('flux_limit.data', output_dir + '/flux_limit.data')
SU.copy('limiter_drawing.data', output_dir + '/limiter_drawing.data')
SU.copy('fort.8', output_dir + '/fort.8')


os.remove('./flux.data')
os.remove('./flux_limit.data')
os.remove('./limiter_drawing.data')
os.remove('./fort.8')

#
