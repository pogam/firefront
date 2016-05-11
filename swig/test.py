import sys
sys.path.append('../')
import forefire, os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import pdb 
import imp 
from scipy.io import netcdf


###################################################
def get_dimension_name_and_value(key,nc_variable):
    
    dim_var  = np.array(                                 nc_variable[key].keys()  )
    dim_var_ = np.array([dimension[:-1] for dimension in nc_variable[key] ])
    idx_dim_var = np.where(dim_var_ == 'dim')

    dimensions_dimX = tuple(sorted(dim_var[idx_dim_var])[::-1])
    dimensions_name = []
    for dim in dimensions_dimX:
        dimensions_name.append(nc_variable[key][dim])

    nco.createVariable(key,nc_variable[key]['type'], dimensions_name )

    dimensions_value = [ nco.dimensions[dim] for dim in dimensions_name ]
    
    return  dimensions_name, dimensions_value


###################################################
def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z


#################################
def getLocationFromLine(line):

    llv = line.split("loc=(")
    if len(llv) < 2: 
        return None
    llr = llv[1].split(",");
    if len(llr) < 3: 
        return None
    return (float(llr[0]),float(llr[1]))


#################################
def dist(a,b):
    return np.sqrt(np.power(a[0]-b[0],2)+np.power(a[1]-b[1],2))


#################################
def printToPathe(linePrinted):
    fronts = linePrinted.split("FireFront")
    pathes = []
    for front in fronts[1:]:
        nodes =front.split("FireNode")[1:]
        if len(nodes) > 0: 
            Path = mpath.Path
           
            codes = []
            verts = []
            lastNode = getLocationFromLine(nodes[0])
            firstNode = getLocationFromLine(nodes[0])
            
 
            codes.append(Path.MOVETO)
            verts.append(firstNode)      

            for node in nodes[:]:
                newNode = getLocationFromLine(node)
                codes.append(Path.LINETO)
                verts.append(newNode)         
                lastNode = newNode
                
    
            codes.append(Path.LINETO)
            verts.append(firstNode)          
           
            pathes.append(mpath.Path(verts, codes))

    return pathes;


###################################################
###################################################


# dimension
nx = 52
ny = 52
#load attribute
attribute = imp.load_source('na','./parametersProperties.py')

#add domain length to attribute and get atmo reso
Lx = 1.*(attribute.domain['NEx']-attribute.domain['SWx'])
Ly = 1.*(attribute.domain['NEy']-attribute.domain['SWy'])
attribute.domain['Lx'] = Lx
attribute.domain['Ly'] = Ly
atmo_dx = Lx/nx
atmo_dy = Ly/nx
y_center = .5 * (Ly - 2 * atmo_dy) #of the domain without the grid pt on the side

#init param
windU_in = 10 # cte wind speed
windV_in = 0
mydatanc = './test/mydata.nc'
mybmapnc = './test/mybmap.nc'
x_location_fireLine = 100.   # location of a fix fire line defined in the bmap file (m)
width_fireLine      = 200.
burningDuration     = 100.   # (s)

#this need to be defined here for bmap matrix
minimalPropagativeFrontDepth = 0.1
spatialIncrement             = .3
BMapsResolution = max(spatialIncrement/np.sqrt(2), minimalPropagativeFrontDepth)


#---------------------------------------------
#create data nc file for ForeFire Layer Models
#---------------------------------------------
print 'write data file: ', mydatanc.split('/')[-1]
fnameout = mydatanc
nc_dimension = {'DIMX': -999, 'DIMY': -999, 'DIMZ': -999, 'DIMT': -999, 'domdim': -999}


var4dim = {'dim1': 'DIMX',   'dim2': 'DIMY',  'dim3': 'DIMZ', 'dim4': 'DIMT', 'type': 'float' }
nc_variable  = {'vaporFlux' : var4dim, \
                'heatFlux'  : var4dim, \
                'fuel'      : var4dim, \
                'windU'     : var4dim, \
                'windV'     : var4dim, \
                'altitude'  : var4dim, \
                'parameters': {'dim1': 'domdim', 'type': 'c'},  
                'domain'    : {'dim1': 'domdim', 'type': 'c'}  }

nc_attribute  = {'vaporFlux': attribute.vaporFlux,
                'heatFlux'  : attribute.heatFlux, \
                'fuel'      : attribute.fuel, \
                'windU'     : attribute.windU, \
                'windV'     : attribute.windV, \
                'altitude'  : attribute.altitude, \
                'parameters': attribute.parameters,                
                'domain'    : attribute.domain}


#open nc file
nco = netcdf.netcdf_file(fnameout, 'w')

#set up dimension
for dim in nc_dimension.keys():
    if dim == 'DIMX':
        nco.createDimension(dim, nx )
    elif dim == 'DIMY':
        nco.createDimension(dim, ny )
    elif dim == 'DIMZ':
        nco.createDimension(dim, 1 )
    elif dim == 'DIMT':
        nco.createDimension(dim, 1 )
    elif dim == 'domdim':
        nco.createDimension(dim, 1 )
    else:
        print 'shoul not be here'
        pdb.set_trace()

#create variables
for key in nc_variable:
    dimensions_name, dimensions_value = get_dimension_name_and_value(key,nc_variable)
    nco.variables[key][:] = np.zeros(dimensions_value)

    if   key == 'windU': 
        nco.variables[key][:] = np.zeros(dimensions_value) + windU_in
    elif key == 'windV': 
        nco.variables[key][:] = np.zeros(dimensions_value) + windV_in
    elif key == 'vaporFlux': 
        nco.variables[key][:] = np.zeros(dimensions_value) + 3

    else:
        nco.variables[key][:] = np.zeros(dimensions_value)

    for attvar in nc_attribute[key]: 
        nco.variables[key]._attributes[attvar] =  nc_attribute[key][attvar]
         

#close nc file
nco.close()


#---------------------------------------------
# create bmap nc file to force fire propapation
#---------------------------------------------
print 'write bmap file: ', mybmapnc.split('/')[-1]
fnameout = mybmapnc
nc_dimension = {'DIMX': -999, 'DIMY': -999, 'C_DIMX': -999, 'C_DIMY': -999, 'domdim': -999}

var2dim_atm  = {'dim1': 'C_DIMX',   'dim2': 'C_DIMY',  'type': 'int32' }
var2dim_fire = {'dim1': 'DIMX',     'dim2': 'DIMY',    'type': 'float' }

nc_variable  = {'arrival_time_of_front' : var2dim_fire, \
                'cell_active'           : var2dim_atm,  
                'domain'                : {'dim1': 'domdim', 'type': 'c'}  }

nc_attribute  = {'arrival_time_of_front': None,
                'cell_active'           : None, \
                'domain'                : merge_two_dicts(attribute.domain,attribute.parameters)}

#sys.exit()
#open nc file
nco = netcdf.netcdf_file(fnameout, 'w')

#set up dimension
for dim in nc_dimension.keys():
    if dim == 'C_DIMX':
        nco.createDimension(dim, nx )
    elif dim == 'C_DIMY':
        nco.createDimension(dim, ny )
    elif dim == 'DIMX':
        nbre_firePt_per_atmoCell_x = (int(atmo_dx/BMapsResolution) + 1 ) 
        nco.createDimension(dim, nbre_firePt_per_atmoCell_x * (nx))
        bmpa_dx = atmo_dx / nbre_firePt_per_atmoCell_x
        grid_x = np.arange(attribute.domain['SWx'],attribute.domain['NEx']+bmpa_dx, bmpa_dx)
    elif dim == 'DIMY':
        nbre_firePt_per_atmoCell_y = (int(atmo_dy/BMapsResolution) + 1 ) 
        nco.createDimension(dim, nbre_firePt_per_atmoCell_y * (ny))
        bmpa_dy = atmo_dy / nbre_firePt_per_atmoCell_y
        grid_y = np.arange(attribute.domain['SWy'],attribute.domain['NEy']+bmpa_dy, bmpa_dy)
    elif dim == 'domdim':
        nco.createDimension(dim, 1 )
    else:
        print 'shoul not be here'
        pdb.set_trace()

#create variables
for key in nc_variable:
    dimensions_name, dimensions_value = get_dimension_name_and_value(key,nc_variable)
    nco.variables[key][:] = np.zeros(dimensions_value)

    if   key == 'arrival_time_of_front': 
            destAT = nco.variables['arrival_time_of_front']
            x_line =  x_location_fireLine
            y_line =  width_fireLine
            i_x_line = (np.abs(grid_x-x_line)).argmin()
            j_y_line = np.where( (grid_y>=(y_center-.5*y_line)) & (grid_y<=(y_center+.5*y_line)) )
            idx_line_start = (np.zeros_like(j_y_line[0])+i_x_line, j_y_line[0] )
            values = np.zeros(dimensions_value[::-1]) -9999
            values[idx_line_start] = 0
            destAT[:] = values.T

    elif key == 'cell_active': 
        nco.variables[key][:] = np.ones(dimensions_value,dtype=np.int) 
    
    else:
        nco.variables[key][:] = np.zeros(dimensions_value)

    if nc_attribute[key] is not None:
        for attvar in nc_attribute[key]: 
            nco.variables[key]._attributes[attvar] =  nc_attribute[key][attvar]
         

#close nc file
nco.close()


#---------------------------------------------
# set ForeFire simulation
#---------------------------------------------
ff = forefire.PLibForeFire()

#set param
ff.setString('ForeFireDataDirectory','test')
ff.setString('fireOutputDirectory','Output')
ff.setInt('outputsUpdate',2)
ff.setString('NetCDFfile',    mydatanc.split('/')[-1])
ff.setString('fluxNetCDFfile',mydatanc.split('/')[-1])
ff.setString('fuelsTableFile','fuels.ff')
ff.setString('BMapFiles', mybmapnc.split('/')[-1])
ff.setDouble("spatialIncrement",spatialIncrement)
ff.setDouble("perimeterResolution",1)
ff.setDouble("minimalPropagativeFrontDepth",minimalPropagativeFrontDepth)
ff.setDouble("nominalHeatFlux",100000)
ff.setDouble("nominalVaporFlux",0.2)
ff.setDouble("burningDuration",burningDuration)

ff.setDouble("bmapOutputUpdate",1)


#set domain
ff.setInt("atmoNX",nx)
ff.setInt("atmoNY",ny)
ff.execute("FireDomain[sw=(%f,%f,0.);ne=(%f,%f,0.);t=0.]"%(attribute.domain['SWx'],attribute.domain['SWy'],attribute.domain['NEx'],attribute.domain['NEy']))


#set propagation model
ff.addLayer("propagation","TroisPourcent","propagationModel")

print "resolution of bmap is ", ff.getString("bmapResolution")

#set fire line
ff.execute("\tFireFront[t=0.]")
ff.execute("\t\tFireNode[loc=(500,600,0.);vel=(0.,0.,0.);t=0.]")
ff.execute("\t\tFireNode[loc=(500,400,0.);vel=(0.,0.,0.);t=0.]")
ff.execute("\t\tFireNode[loc=(490,400,0.);vel=(0.,0.,0.);t=0.]")
ff.execute("\t\tFireNode[loc=(490,600,0.);vel=(0.,0.,0.);t=0.]")


#---------------------------------------------
# run ForeFire simulation
#---------------------------------------------
pathes = []
step = 1
N_step = 10
for i in np.arange(1,N_step):
    print "goTo[t=%f]"%(i*step)
    ff.execute("goTo[t=%f]"%(i*step))
    pathes += printToPathe( ff.execute("print[]"))


#---------------------------------------------
# plot ForeFire perimeter
#---------------------------------------------
fig, ax = plt.subplots()
 
#tab = np.transpose(ff.getDoubleArray("BMap"))[0]
#CS = ax.imshow(tab, origin='lower', cmap=plt.cm.gray, interpolation='nearest',\
#               extent=(0,tab.shape[1]*np.float(ff.getString("bmapResolution")),0,tab.shape[0]*np.float(ff.getString("bmapResolution"))))
#cbar = plt.colorbar(CS)
#cbar.ax.set_ylabel('v')
cmap = mpl.cm.get_cmap('jet')
for i, path in enumerate(pathes):
    rgba = cmap(1.*i/(N_step-1))
    patch = mpatches.PathPatch(path,edgecolor=rgba, facecolor='none', alpha=1)
    ax.add_patch(patch)

ax.grid()
ax.axis('equal')
ax.set_xlim(attribute.domain['SWx']+atmo_dx,attribute.domain['NEx']-atmo_dx)
ax.set_ylim(attribute.domain['SWy']+atmo_dy,attribute.domain['NEy']-atmo_dy)
plt.show()


