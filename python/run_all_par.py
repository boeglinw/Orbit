##################Importing Modules##################
import os
from LT import parameterfile as PF
import LT.box as B
import numpy as np

# parse arguments

import argparse as AG

##################Defining Constants################

# Convert mm to m

mm2m = 1.e-3

# Convert degrees to radians

dtr = np.pi/180

####################Defining Functions##############
#rotate positions of detectors
def posi_rot(x,y,phi): #rotation around r
    newX=x*np.cos(phi)+y*np.sin(phi)
    newY=-x*np.sin(phi)+y*np.cos(phi)
    return np.array([newX,newY])
# get a polar angle for a 2d vector
def pol_angle(rx,ry):
    cphi = rx/np.sqrt(rx**2+ry**2)
    phi = np.arccos(cphi)
    phic = np.where(ry>0, phi, 2.*np.pi - phi) # set all values where rx  < 0 to 2pi - phi
    return phic

    # return the angle corresponding to these
def rotate_det(phi, theta, alpha):
    import pdb
    eps = 1.e-7
    phi_r = np.zeros_like(phi)
    theta_r = np.zeros_like(theta)
    # unit vector of detector orientation
    nx = -np.cos(phi)*np.sin(theta)
    ny = -np.sin(phi)
    nz = np.cos(phi)*np.cos(theta)
    n = np.array([nx,ny,nz])
    # rotated vector
    nxr = nx
    nyr = ny*np.cos(alpha) - nz*np.sin(alpha)
    nzr = ny*np.sin(alpha) + nz*np.cos(alpha)
    nr = np.array([nxr, nyr, nzr])
    # now claculate new angles
    # handle limits
    case1 = (np.abs(nzr) <= eps)&(nyr>0.)
    case2 = (np.abs(nzr) <= eps)&(nyr <0.)
    case3 = (np.abs(nzr) <= eps)&( np.abs(nyr) <= eps)&(nx > 0.)
    case4 = (np.abs(nzr) <= eps)&( np.abs(nyr) <= eps)&(nx < 0.)
    case5 = np.logical_not ( case1 | case2 | case3 | case4)
    #pdb.set_trace()
    # if (np.abs(nzr) <= eps)&(nyr>0.):
    theta_r[case1] = np.pi/2.
    phi_r[case1] = 3./2.*np.pi + theta[case1]
    # elif(np.abs(nzr) <= eps)&(nyr <0.):
    theta_r[case2] =  np.pi/2.
    phi_r = np.pi/2. - theta
    # elif(np.abs(nzr) <= eps)&( np.abs(nyr) <= eps):
    phi_r[case3]  = 0.
    # if (nx > 0.):
    theta_r[case3] = 3.*np.pi/2.
    #    else:
    phi_r[case4] = 0.
    theta_r[case4] = np.pi/2.
    #else:
    theta_r[case5] = -np.arctan(nxr[case5]/nzr[case5])
    sphi = -nyr[case5]
    cphi = nzr[case5]/np.cos(theta_r[case5])
    phi_r[case5] = pol_angle(cphi, sphi)
    return np.array([phi_r, theta_r])
    
def phi_change(phi,theta,alpha):
    A = np.sin(phi)*np.cos(alpha)
    B = np.cos(phi)*np.cos(theta)*np.sin(alpha)
    return np.arcsin(A-B)/dtr

def theta_change(phi,theta,alpha):
    x_prime = -1.0*np.cos(phi)*np.sin(theta)
    z_prime = np.sin(phi)*np.sin(alpha)+np.cos(phi)*np.cos(theta)*np.cos(alpha)
    A = z_prime/(np.sqrt(x_prime**2 + z_prime**2))
    return np.arccos(A)/dtr

#########################Proper Files###############

# reverse direction
reverse = False

parser = AG.ArgumentParser()
parser.add_argument("control_file", nargs = '?', help="Control file ", default = 'control.data')
args = parser.parse_args()

# open control file
c_file = args.control_file

cd = PF.pfile(c_file)

machine = cd.get_value('machine')

input_root = cd.get_value('input_root')
input_dir_ext = cd.get_value('input_dir_ext')

nml_dir = input_root + machine + input_dir_ext + '//'

#nml_filename = 'MAST_p6'
nml_filename = cd.get_value('input_file')

make_new_nml = cd.get_value('make_new_nml', var_type = cd.Bool)


#static_dir = 'MAST_input'
static_dir = cd.get_value('static_dir')

#dynamic_dir = 'MAST_input'
dynamic_dir = cd.get_value('dynamic_dir')

dyn_file = cd.get_value('dynamic_file')
stat_file = cd.get_value('static_file')

# see of one wants reverse the particle direction (usefule for determining detector locations)
try:
    reverse = cd.get_value('reverse_direction', var_type = cd.Bool)
except:
    print 'no reverse_direction found use default value : False'
    reverse = False

if reverse:
    print 'Changed particle direction, IFOR  in the static file need to be set to 1!'
else:
    print 'Normal particle direction, IFOR  in the static file need to be set to 0!'
    
# Total number of detectors in the system

det = cd.get_value('N_det', var_type = int)

nml_ext = cd.get_value('nml_ext')

orbit_command = cd.get_value('orbit_command') + ' ' + c_file

plot_command = cd.get_value('plot_command') + ' ' + c_file

##########################Obtaining Parameters and Data############################


# Grabbing dynamic data file

# dd = B.get_file('./'+ dynamic_dir + '/dynamic_input.data') 
dd = B.get_file(dynamic_dir + '/' + dyn_file) 

dpar = dd.par        
                        
# Number of detectors used in the analysis

detectors = dd.par.get_value('detectors',var_type = int)

# Which detectors are used in the calculations

dn = np.array(dpar.get_value('detector_number').split(','), dtype = int)

# Port angle of each detector when the Reciprocating probe arm is at a rotation angle of 0

# Port angle was input in degrees

ppb = B.get_data(dd, 'phi_port_base')*dtr
if reverse:
    # to change direction of particle add 180 degrees
    ppb += np.pi

# Theta angle of each detector when the Reciprocating probe arm is at a rotation angle of 0

# Port angle was input in degrees

tpb = B.get_data(dd,'theta_port_base')*dtr

# Horizontal offset, from the center of the diagnostic, of each detector 

# Horizontal offset was input in mm

dxo = B.get_data(dd, 'detector_horizontal_offset')*mm2m

# Radial offset of each detector (measured from the base)

# Radial offset was input in mm

dro = B.get_data(dd, 'detector_radial_axis_offset')*mm2m

# Height offset of each detector (measured from the center)

# Height offset was input in mm

dyo = B.get_data(dd, 'detector_height_offest')*mm2m

# get the assigned channel numbers
channel_number = B.get_data(dd, 'ch')

# get the assigned channel numbers
detector_id = B.get_data(dd, 'detector_id')

# Alpha is the rotational angle of the RP arm 

# Alpha was input in degrees

alpha = dpar.get_value('RP_rotation')*dtr

####################Making Calculations##########################

# ZDist is the height of each detector relative to the center of the tokamak

ZDist = dpar.get_value('ZDist') + posi_rot(dxo,dyo,alpha)[1]


# RDist is the distance from the center of the tokamak of each detector in meters

RDist =  dpar.get_value('RDist') + dpar.get_value('RDist_offset')+ dro



# phdangle is the toroidal angle of each detector relative to the 0 angle
# phangle is in radians

phdangle = dpar.get_value('PHDangle')*dtr + np.arctan(posi_rot(dxo,dyo,alpha)[0]/RDist) 
  
# The following are the calculations for the change in phi and theta angles
# for a rotation of the reciprocating probe arm 
        
if alpha == 0.:
    
    phi_port = ppb/dtr
    theta_port = tpb/dtr
    
else:
    phi_port,theta_port = rotate_det(ppb,tpb,alpha)
    phi_port /= dtr
    theta_port /= dtr
        
#    theta_port = theta_change(ppb,tpb,alpha)
#    
#    phi_port = phi_change(ppb,tpb,alpha)
#    
#    theta_port = theta_change(ppb,tpb,alpha)

for i in range(det):
    print 'Phi_port was ',ppb[i]/dtr,' degrees. After a rotation of ', alpha/dtr,'degrees, it is now ',phi_port[i]
    print 'Theta was ', tpb[i]/dtr,'degrees. After a rotation of ',alpha/dtr,' degrees, it is now ', theta_port[i]

######################Creating New Dictionary########################

# Physical Parameter Dictionary

PA = {\
'PHDangle':phdangle,\
'RDist':RDist,\
'phi_port': phi_port,\
'theta_port':theta_port,\
'channel_number':channel_number,\
'detector_id':detector_id,\
'ZDist': ZDist,\
'detector_number':dn,\
'detectors':detectors\
}
    
PAK = PA.keys()
PAG = PA.get

######################Make a new NML file#########################
if make_new_nml == True:
    
    staticf =  open(static_dir + '/' + stat_file, 'r').readlines() #opens and reads static file
    
    bothf = open(nml_dir + nml_filename + nml_ext,'w+')   #creates nml writable file 
    
    # Writing New nml file    
    
    bothf.writelines(staticf[:staticf.index(' &orbit_par\n')+2])     #writes the static file into the new nml file     
    
    # Selecting the inputs under orbitpar
    
    orbit_par = ['detector_number', 'channel_number', 'detector_id', 'detectors', 'theta_port', 'phi_port']
    
    # Writing Parameters into orbit_par in nml file
    
    for i in orbit_par:
        bothf.write('\n\n')
        for k in range(det):
            try:
                bothf.write('  ' + i + '({0:1})'.format(k+1) + '= {0:1}\n'.format(PAG(i)[k]))
            except:
                bothf.write('  ' + i + '= {0:1}\n'.format(PAG(i)))
                break        
    bothf.writelines(staticf[staticf.index(' &orbit_par\n')+1:staticf.index(' &detection\n')+2])               

    # Choosing keys for detection
    
    detection = []
    for i in PAK:
        if i not in orbit_par:
            detection.append(i)
    
    # Writing Parameters into detection in nml file
    
    for i in detection:
        bothf.write('\n\n')
        for k in range(det):
            try:
                bothf.write('  ' + i + '({0:1})'.format(k+1) + '= {0:1}\n'.format(PAG(i)[k]))
            except:
                bothf.write('  ' + i + '= {0:1}\n'.format(PAG(i)))
                break
    
    bothf.writelines(staticf[staticf.index(' &detection\n')+1:])              
    bothf.close()



os.system(orbit_command)


# about plotting
answer = False

while answer == False:

	a  = raw_input('Do you want to view the plots? (yes/no)')

	if a == 'yes':

		os.system(plot_command)
		answer = True

	elif a == 'no':

		print '\n Please continue reading below. \n'
		answer = True

	else: 

		print 'Please enter yes or no.'

"""
# Determing what to enter in the RP Remote Control

SpGOn = lambda rp, pdo: 1200 - (2526 - rp*1000. - pdo*1000.)

SpGOff = lambda rp, pdo: 1200 - (2526 - rp*1000. - 105.4 - pdo*1000)

PD_offset = (dro.max() + 42)*1.e-3 - .165

RP_min = PAG('RDist').min()

print '\nWhen the collimator entrance closest to the center of the MAST Tokamak is at', RP_min,' meters:\n'
print '\n- The outer edge of the shell is ',RP_min- 11.e-3,' meters from the center of the plasma.\n'
print '- The value that must be inputted into the remote RP controller can be one of the following two:\n'
print '-- If the gas is on, input ',SpGOn(RP_min,PD_offset),'\n'
print '-- If the gas is off, input ',SpGOff(RP_min,PD_offset),'\n'
print 'Please remember to stay between 248 and 500!\n'
print 'Never input a value above 500!'

print '\nHave a nice day.'
"""
