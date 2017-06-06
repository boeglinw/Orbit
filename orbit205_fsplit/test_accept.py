import numpy as np

xc, xd, d, xal = input('enter xc, xd, d, xal :')
print 'entered : ', xc, xd, d, xal
coll_size = 2.*xc
det_size = 2.*xd
h_offset = xal
Dist = d

if (coll_size <= det_size):
    full_opening = coll_size
    theta1n = np.arctan( (coll_size - det_size - 2.*h_offset)/(2.*Dist) )
    theta1p = np.arctan( (det_size - coll_size - 2.*h_offset)/(2.*Dist) )
else:
    full_opening = det_size
    theta1p = np.arctan( (coll_size - det_size - 2.*h_offset)/(2.*Dist) )
    theta1n = np.arctan( (det_size - coll_size - 2.*h_offset)/(2.*Dist) )
theta2n = np.arctan( -(det_size + coll_size + 2.*h_offset)/(2.*Dist) )
theta2p = np.arctan( (det_size + coll_size - 2.*h_offset)/(2.*Dist) )
# acceptance
acc_1 = full_opening * ( np.sin(theta1p) - np.sin(theta1n) ) 
acc_2 =     ((det_size + coll_size + 2.*h_offset)/2.) * (np.sin(theta1n) - np.sin(theta2n) )
acc_3 =     -Dist*(np.cos(theta1n) - np.cos(theta2n))
acc_4 =     ((det_size + coll_size - 2.*h_offset)/2.) * (np.sin(theta2p) - np.sin(theta1p) )
acc_5 =     Dist*(np.cos(theta2p) - np.cos(theta1p))

acc_f = full_opening * ( np.sin(theta1p) - np.sin(theta1n) ) + \
    ((det_size + coll_size + 2.*h_offset)/2.) * (np.sin(theta1n) - np.sin(theta2n) ) - \
    Dist*(np.cos(theta1n) - np.cos(theta2n)) +\
    ((det_size + coll_size - 2.*h_offset)/2.) * (np.sin(theta2p) - np.sin(theta1p) ) + \
    Dist*(np.cos(theta2p) - np.cos(theta1p))

print "acc_1 = ",acc_1
print "acc_2 = ", acc_2
print "acc_3 = ", acc_3
print "acc_4 = ", acc_4
print "acc_5 = ", acc_5
print 'acceptance', acc_f
print 'sum(acc_i) = ', acc_1 + acc_2 + acc_3 + acc_4 + acc_5
