# Importing Modules
import run_command as rc
import shutil as SU
import os
import sys
import LT.box as B
import numpy as np

# Proper Files

input_file = 'MAST_input'

input_ext = '.nml'

input_name = input_file +input_ext

# Obtaining Parameters and Data


dynamic_data = B.get_file('../'+ input_file + '/dynamic_input.data') 

dynamic_parameters = dynamic_data.par.get_all_data()

for k in dynamic_parameters.keys():
    try:
        dynamic_parameters[k]=float(dynamic_parameters[k])
    except:
        print "\n Cannot make a float  out of this : ", dynamic_parameters[k], k,'\n'
        
detectors = 4 # Total number of detectors in the system

try:        
	dynamic_parameters['detector_number']=np.array(dynamic_parameters['detector_number'].split(','), dtype = int)
except:
	dynamic_parameters['detector_number'] = [int(dynamic_parameters['detector_number'])]

dn = dynamic_parameters['detector_number']

# location of the base

phi_port_base = B.get_data(dynamic_data, 'phi_port_base')

# position offsets for the individual detectors relative to the base

detector_horizontal_offset = B.get_data(dynamic_data, 'detector_horizontal_offset')

detector_radial_axis_offset = B.get_data(dynamic_data, 'detector_radial_axis_offset')

detector_height_offset = B.get_data(dynamic_data, 'detector_height_offest')

# Making Calculations

dynamic_parameters['ZDist'] = dynamic_parameters['ZDist'] + detector_height_offset*1.e-3

# With respect to the radial position of the RP

dp =  []

for i,m in enumerate(detector_radial_axis_offset):
	dp.append(detector_radial_axis_offset.max()-detector_radial_axis_offset[i])

dp = np.array(dp)

dynamic_parameters['RDist'] =  np.round(dynamic_parameters['RDist'] + dp*1.e-3,2)

phdangle = []

for i in range(detectors):
    phdangle.append(np.round(dynamic_parameters['PHDangle'] + np.arctan(detector_horizontal_offset[i]*1.e-3/dynamic_parameters['RDist'][i]),2))

dynamic_parameters['PHDangle'] = np.array(phdangle) 
  
# With respect to the angular change of the RP  

if dynamic_parameters['RP_rotation'] == 0.:
    
    phi_port = phi_port_base
    dynamic_parameters['theta_port'] = [dynamic_parameters['theta_port']]*detectors
    
else:
    
    phi_port = []
    theta_port = []
    
    for i in range(detectors):
        
        # Phi, Theta, and RP_rotation angles.
        
        phi = np.radians(phi_port_base[i])
        alpha = np.radians(dynamic_parameters['RP_rotation'])
        theta = np.radians(dynamic_parameters['theta_port'])
        
        # New Phi port Calculation
               
        phi_port.append(np.degrees(np.arcsin(np.sin(phi)*np.cos(alpha)-np.cos(phi)*np.cos(theta)*np.sin(alpha))))
        print 'Phi_port was ',np.degrees(phi),' degrees. After a rotation of ', np.degrees(alpha),', it is now ',phi_port[i]
        
        # New Theta port Calculation
        
        x_prime = -1.0*np.cos(phi)*np.sin(theta)
        z_prime = np.sin(phi)*np.sin(alpha)+np.cos(phi)*np.cos(theta)*np.cos(alpha)
        theta_port.append(np.degrees(np.arccos(z_prime/(np.sqrt(x_prime**2 + z_prime**2)))))
    
        print 'Theta was ', np.degrees(theta),'degrees. After a rotation of ',np.degrees(alpha),' degrees, it is now ', theta_port[i]
    
    phi_port = np.array(phi_port)
    dynamic_parameters['theta_port'] = np.array(theta_port)
    
    
# Opening the NML files

staticf =  open('./'+ input_file + '/Static_MAST.nml', 'r').readlines() #opens and reads static file

os.chdir('./'+input_file)

bothf = open('MAST_p6.nml','w+')             #creates MAST_p6.nml writable file (will eventually become Mast_p5.nml)

   
# Writing New nml file    

bothf.writelines(staticf[:staticf.index(' &orbit_par\n')+2])     #writes the static file into the new nml file     

orbit_par = [dynamic_parameters.keys()[i] for i in [dynamic_parameters.keys().index('detectors'),dynamic_parameters.keys().index('detector_number'),dynamic_parameters.keys().index('theta_port')]];orbit_par.append('phi_port')

# print orbit_par

#bothf.write('\n &orbit_par \n\n')

for i in orbit_par:
    if i == 'phi_port':
        orbit_par.remove(i)
        for k in range(detectors):
            if k == len(range(detectors))-1:
                bothf.write('  ' + i + '({0:1})'.format(k+1) + '= {0:1}\n\n'.format(phi_port[dn[k]-1])) 
            else:  
                bothf.write('  ' + i + '({0:1})'.format(k+1) + '= {0:1}\n'.format(phi_port[dn[k]-1]))     
    
    elif i == 'detectors':
        bothf.write('  ' + i + '= {0:1}\n\n'.format(dynamic_parameters[i]))
			
    
    else:
        for k in range(detectors):
                if i == 'detector_number':
                    if k == len(range(detectors))-1:
                        bothf.write('  ' + i +'({0:1})'.format(k+1) + '= {0:1}\n\n'.format(dynamic_parameters[i][k]))
                    else:
                        bothf.write('  ' + i + '({0:1})'.format(k+1) + '= {0:1}\n'.format(dynamic_parameters[i][k]))        
                elif k == len(range(detectors))-1:
                	   bothf.write('  ' + i +'({0:1})'.format(k+1) + '= {0:1}\n\n'.format(dynamic_parameters[i][dn[k]-1]))
            	else:
                   bothf.write('  ' + i + '({0:1})'.format(k+1) + '= {0:1}\n'.format(dynamic_parameters[i][dn[k]-1]))        
                    
bothf.writelines(staticf[staticf.index(' &orbit_par\n')+1:staticf.index(' &detection\n')+2])               

detection = dynamic_parameters.keys(); detection.remove('RP_rotation')
for i in orbit_par:
    try:
        detection.remove(i)
    except:
        print i,'could not be removed.'    

# print detection

for i in detection:
        for k in range(detectors):
            if k == len(range(detectors))-1:
                bothf.write('  ' + i +'({0:1})'.format(k+1) + '= {0:1}\n\n'.format(dynamic_parameters[i][dn[k]-1]))
            else:
                bothf.write('  ' + i + '({0:1})'.format(k+1) + '= {0:1}\n'.format(dynamic_parameters[i][dn[k]-1]))
                
bothf.writelines(staticf[staticf.index(' &detection\n')+1:])              
bothf.close()

os.chdir('..')

os.system('python ./python/run_orbit_nml.py')

answer = False

while answer == False:

	a  = raw_input('Do you want to view the plots? (yes/no)')

	if a == 'yes':

		os.system('python ./python/plot_orbits_combined.py')
		answer = True

	elif a == 'no':

		print '\n Please continue reading below. \n'
		answer = True

	else: 

		print 'Please enter yes or no.'

# Determing what to enter in the RP Remote Control

SpGOn = lambda rp, pdo: 1200 - (2526 - rp*1000. - pdo*1000.)

SpGOff = lambda rp, pdo: 1200 - (2526 - rp*1000. - 105.4 - pdo*1000)

PD_offset = (detector_radial_axis_offset.max() + 42)*1.e-3 - .165

RP_min = dynamic_parameters['RDist'].min()

print '\nWhen the collimator entrance closest to the center of the MAST Tokamak is at', RP_min,' meters:\n'
print '\n- The outer edge of the shell is ',RP_min- 11.e-3,' meters from the center of the plasma.\n'
print '- The value that must be inputted into the remote RP controller can be one of the following two:\n'
print '-- If the gas is on, input ',SpGOn(RP_min,PD_offset),'\n'
print '-- If the gas is off, input ',SpGOff(RP_min,PD_offset),'\n'
print 'Please remember to stay between 248 and 500!\n'
print 'Never input a value above 500!'

#if SpGOn(RP_min,PD_offset) >= 500. or SpGOff(RP_min,PD_offset) >=500.:
#	if SpGOn(RP_min,PD_offset) >= 500.:
#		print '\nWARNING THE SET POINT WITH THE GAS ON IS AT OR ABOVE 500!'
#	elif SpGOff(RP_min,PD_offset) >= 500.:
#		print '\nWARNING THE SET POINT WITH THE GAS OFF IS AT OR ABOVE 500!'

print '\nHave a nice day.'



	
