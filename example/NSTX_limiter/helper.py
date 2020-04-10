# prepare geometry data files

import numpy as np
import LT.box as B


# fill these values interactively
n = 0

comments = []
phi = []
rs = []
zs = []

rl = []
zl = []

# inner side, same for all toroidal angle zones

rs = np.array([ 0.61706, 0.57150, 0.57150, 0.27940, 0.27940, 0.18515, 
  0.18515, 0.27940, 0.27940, 0.57150, 0.57150, 0.61706])

zs = np.array([-1.62801,-1.62800,-1.60338,-1.60338,-1.17137,-1.00812, 
  1.00812, 1.17137, 1.60338, 1.60338, 1.62800, 1.62801])


header_comment = []

header_comment.append('! number of toroidal regions')
header_comment.append('NSTX passive stabilizer plates, divertor plates, and RF antenna. starting angle of first toroidal region (L-I bays)  ')


# outer zone
phi.append(0.)
comments.append('! NB armor (mid I-H bays) per FDR data 1NOV99 (15FEB00)')
rl.append(np.array([0.61706, 1.04327, 
  1.04327, 1.31915, 1.33584, 1.33584, 1.48509, 1.69690, 
  1.69690, 1.48509, 1.33584, 1.33584, 1.31915, 1.04327, 
  1.04327, 0.61706]))
zl.append(np.array([-1.62801,-1.46012, 
 -1.43003,-1.03967,-1.03966,-0.99761,-0.54501,-0.54450, 
  0.54450, 0.54501, 0.99761, 1.03966, 1.03967, 1.43003, 
  1.46012, 1.62801]))


phi.append(8.315)
comments.append('! Proton detector port')
rl.append(np.array([0.61706, 1.04327, 
  1.04327, 1.31915, 1.33584, 1.33584, 1.48509, 1.69690, 
  1.69690, 1.8,     1.8,     1.69690, 1.69690,
  1.48509, 1.33584, 1.33584, 1.31915, 1.04327, 
  1.04327, 0.61706]))
zl.append(np.array([-1.62801,-1.46012, 
 -1.43003,-1.03967,-1.03966,-0.99761,-0.54501,-0.54450, 
  0.20,    0.201,    0.3,     0.301,     0.54450, 
  0.54501, 0.99761, 1.03966, 1.03967, 1.43003, 
  1.46012, 1.62801]))

phi.append(11.685)
comments.append('! continuation of first toroidal region (L-I bays)')
rl.append(np.array([0.61706, 1.04327, 
  1.04327, 1.31915, 1.33584, 1.33584, 1.48509, 1.69690, 
  1.69690, 1.48509, 1.33584, 1.33584, 1.31915, 1.04327, 
  1.04327, 0.61706]))
zl.append(np.array([-1.62801,-1.46012, 
 -1.43003,-1.03967,-1.03966,-0.99761,-0.54501,-0.54450, 
  0.54450, 0.54501, 0.99761, 1.03966, 1.03967, 1.43003, 
  1.46012, 1.62801]))
  
  
phi.append(105.0)
comments.append('! NB armor (mid I-H bays) per FDR data 1NOV99 (15FEB00)')
rl.append(np.array([0.61706, 1.04327, 
  1.04327, 1.31915, 1.33584, 1.33584, 1.48509, 1.58433, 
  1.58433, 1.48509, 1.33584, 1.33584, 1.31915, 1.04327, 
  1.04327, 0.61706]))
zl.append(np.array([-1.62801,-1.46012, 
 -1.43003,-1.03967,-1.03966,-0.99761,-0.54501,-0.54450, 
  0.54450, 0.54501, 0.99761, 1.03966, 1.03967, 1.43003, 
  1.46012, 1.62801]))

phi.append(147.0)
comments.append('! starting angle for H-F bays')
rl.append(np.array([0.61706, 1.04327, 
  1.04327, 1.31915, 1.33584, 1.33584, 1.48509, 1.69690, 
  1.69690, 1.48509, 1.33584, 1.33584, 1.31915, 1.04327, 
  1.04327, 0.61706]))
zl.append(np.array([-1.62801,-1.46012, 
 -1.43003,-1.03967,-1.03966,-0.99761,-0.54501,-0.54450, 
  0.54450, 0.54501, 0.99761, 1.03966, 1.03967, 1.43003, 
  1.46012, 1.62801]))

phi.append(210.0 )
comments.append('! starting angle of HHFW antenna')
rl.append(np.array([0.61706, 1.04327, 
  1.04327, 1.31915, 1.33584, 1.33584, 1.48509, 
  1.4851,  1.5708,  1.5788,  
  1.5820,  1.5788,  1.5708,  1.4851,
  1.48509, 1.33584, 1.33584, 1.31915, 1.04327, 
  1.04327, 0.61706]))
zl.append(np.array([-1.62801,-1.46012, 
 -1.43003,-1.03967,-1.03966,-0.99761,-0.54501,
 -0.5450, -0.1139, -0.0610,  
  0.0000,  0.0610,  0.1139,  0.5450,
  0.54501, 0.99761, 1.03966, 1.03967, 1.43003, 
  1.46012, 1.62801]))

phi.append(300.0)
comments.append('! end of antenna back to bay L')
rl.append(np.array([0.61706, 1.04327, 
  1.04327, 1.31915, 1.33584, 1.33584, 1.48509, 1.69690, 
  1.69690, 1.48509, 1.33584, 1.33584, 1.31915, 1.04327, 
  1.04327, 0.61706]))
zl.append(np.array([-1.62801,-1.46012, 
 -1.43003,-1.03967,-1.03966,-0.99761,-0.54501,-0.54450, 
  0.54450, 0.54501, 0.99761, 1.03966, 1.03967, 1.43003, 
  1.46012, 1.62801]))

o_s = open('small_R.data','w')
header = '# smal R data , common for all zones\n#\n'
# all zones have the same inner structure
header +=  '#! r[f,0]/ z[f,1]/\n'
o_s.write(header)
for i,rv in enumerate(rs):
    o_s.write( '{} {} \n'.format(rv, zs[i]))
o_s.close()

for i,ph in enumerate(phi):
    o_l = open('large_R_{}.data'.format(i),'w')
    header = '#\n'
    header += '#\ phi = {}\n'.format(ph)
    header += '#\ comment = {}\n'.format(comments[i])
    header += '#\n'
    header +=  '#! r[f,0]/ z[f,1]/\n'
    
    o_l.write(header)
    for j, rv in enumerate(rl[i]):
        zv = zl[i][j]
        o_l.write( '{} {} \n'.format(rv, zv))
    o_l.close()
print 'all_done'
