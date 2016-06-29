import numpy as np 
from scipy import interpolate
import sys
import os

###################
def wind_profile(z,u0):
    return u0 * (z/2)**(1./7)

############################
if __name__ == '__main__':
############################

    if not os.path.isfile('./input_cst.py'):
        print 'missing input cstn info file: ', './input_cst.py'
        sys.exit()
    else:
        from input_cst import *

    g = 9.8      # m/s2
    cp = 1005.   # J/kg/K
    R  = 287.058 # J/kg/K
    p00 = 100000.

    #those quantity are defined in  input_cst
    #---------------    
    #date='2000 01 01 0.'
    #nbre_layer    = 10
    #grd_theta     = 285.
    #grd_pressure  = 100000.
    #grd_atl       = 0.
    #top_alt       = 3000.
    #lapse_rate_in = [0.     , 0.     ]
    #u0  = 5.

    level_in      = [grd_atl, top_alt]

    f_lapse_rate = interpolate.interp1d(level_in,lapse_rate_in)
    #grd_T        = grd_theta * (grd_pressure/p00)**(R/cp)
    grd_theta    = grd_T  /  (grd_pressure/p00)**(R/cp)
    
    dry_lapse_rate = g/cp  # K/m  T' = T0 - Gamma_d z'
                           #      T  = T0 - Gamma   z'
                           #      N2 = g/T [Gamma_d - Gamma] 

    prof_alt   = np.linspace(grd_atl, top_alt, nbre_layer)
    lapse_rate = f_lapse_rate(prof_alt[:-1]) 

    prof_u   = wind_profile(prof_alt,u0)
    prof_v   = np.zeros(nbre_layer)
    prof_rh  = np.zeros(nbre_layer)
    prof_T   = grd_T + lapse_rate * (prof_alt[:-1] - grd_atl)
    prof_N   = np.sqrt(g/prof_T * (dry_lapse_rate - np.array(lapse_rate)))

#dump CSTN format for PRE_IDEA.name
    f = open('./CSTN.nam','w')
#date
    f.write(date+'\n')
#nbre level
    f.write('{:d}\n'.format(nbre_layer))
#ground Theta
    f.write('{:3.6f}\n'.format(grd_theta))
#ground pressure
    f.write('{:12.6f}\n'.format(grd_pressure))
#level alt
    f.write( ' '.join(["%.3f" % val for val in prof_alt]) + '\n' )
#level u
    f.write( ' '.join(["%.3f" % val for val in prof_u]) + '\n' )
#level v
    f.write( ' '.join(["%.3f" % val for val in prof_v]) + '\n' )
#level rh
    f.write( ' '.join(["%.3f" % val for val in prof_rh]) + '\n' )
#level N
    f.write( ' '.join(["%.3f" % val for val in prof_N]) + '\n' )
#close file
    f.close()
