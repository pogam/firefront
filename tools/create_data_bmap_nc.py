import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import pdb 
import imp 
from scipy.io import netcdf
import glob
import netCDF4
import parametersProperties as attribute
import f90nml
import datetime
import socket 

if socket.gethostname() == 'pyro':
    path_ForeFire =  '/home/ronan/Src/ForeFire/'
elif socket.gethostname() == 'ubu':
    path_ForeFire =  '/home/paugam/Src/ForeFire/'
else:
    print socket.gethostname() + ' is not defined. stop here'
    sys.exit()

sys.path.append(path_ForeFire)
sys.path.append(path_ForeFire+'/swig/')
sys.path.append(path_ForeFire+'tools/')
try: 
    import forefire
except ImportError:
    print 'could do find ForeFire API'
    print 'check variable path_ForeFire in path_ForeFire/tools/create_data_bmap_nc.py'
    print 'currently, path_ForeFire = ', path_ForeFire
    sys.exit()
##################################################
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


##################################################
def read_MNH_domain( expr_name, dir_mnh = '../02_mnh/'):

    outfiles_mnh = sorted(glob.glob(dir_mnh+expr_name+'*.nc4'))
    if len(outfiles_mnh) == 0:
        print 'missing mnh files'
        pdb.set_trace()

    #read geometry of the MNH domain
    nc = netCDF4.Dataset(outfiles_mnh[1],'r') # get the 001 file

    DeltaY = nc.variables['YHAT'][1]-nc.variables['YHAT'][0]
    DeltaX = nc.variables['XHAT'][1]-nc.variables['XHAT'][0]
    
    MNH_properties= {}
    MNH_properties['nx']  = len(nc.dimensions[u'X'])
    MNH_properties['ny']  = len(nc.dimensions[u'Y'])
    MNH_properties['nz']  = len(nc.dimensions[u'Z'])

    MNH_properties['SWx']  = nc.variables['XHAT'][0]
    MNH_properties['SWy']  = nc.variables['YHAT'][0]
    MNH_properties['SWz']  = 0
    MNH_properties['NEx']  = nc.variables['XHAT'][-1]+DeltaX
    MNH_properties['NEy']  = nc.variables['YHAT'][-1]+DeltaY
    MNH_properties['NEz']  = 0

    MNH_properties['Lx']   = nc.variables['XHAT'][-1]+DeltaX-nc.variables['XHAT'][0]
    MNH_properties['Ly']   = nc.variables['YHAT'][-1]+DeltaY-nc.variables['YHAT'][0]
    MNH_properties['Lz']   = 0
    nc.close()

    nc = netCDF4.Dataset(outfiles_mnh[0],'r')
    MNH_properties['t0']   = nc.variables['DTCUR__TIME'][0]
    MNH_properties['Lt']   = np.Inf

    MNH_properties['date']   = datetime.datetime(*nc.variables['DTCUR__TDATE'][:]).strftime("%Y-%m-%d_%H:%M:%S")

    return MNH_properties


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
def getDomainExtent(line):
    print line
    llv = line.split("sw=(")
    llr = llv[1].split("ne=(");
    return( float( llr[0].split(",")[0]), float(llr[1].split(",")[0]), float(llr[0].split(",")[1]), float(llr[1].split(",")[1]) )


###################################################
if __name__ == '__main__':
###################################################

    ################### ... start input param

    #ff param
    input_ff_param_file = './Inputs/ff-param.nam'

    ################### ... end input param

    #create dir
    ensure_dir('./Inputs/')
    ensure_dir('./ForeFire/')
    ensure_dir('./ForeFire/Outputs/')
    ensure_dir('./MODEL1/')


    #read MNH info
    MNH_namelist = f90nml.read('../02_mnh/EXSEG1.nam')
    expr_name = MNH_namelist['NAM_CONF']['CEXP']
    MNH_domain   = read_MNH_domain(expr_name)

    #read FF info
    nmls = f90nml.read('./Inputs/ff-param.nam')
    
    #ronan ff info
    ronan_param = nmls['FF_RONAN']
    #fix FireLine
    x_location_fireLine = ronan_param['x_location_fireLine']  
    length_fireLine     = ronan_param['length_fireLine'] 
    depth_fireLine      = ronan_param['depth_fireLine'] 
    time_start_fireLine = ronan_param['time_start_fireLine']
    nbre_fire_line      = len(time_start_fireLine)
    #extra perimeter
    flag_run_ff        = ronan_param['flag_run_ff'] 
    flag_extra_contour = ronan_param['flag_extra_contour'] 

    #general ff info
    ff_param = nmls['FF_PARAM']
    burningDuration              = ff_param['burningDuration']
    minimalPropagativeFrontDepth = ff_param['minimalPropagativeFrontDepth']
    spatialIncrement             = ff_param['spatialIncrement']
    perimeterResolution          = ff_param['perimeterResolution']
    nominalHeatFlux              = ff_param['nominalHeatFlux']
    nominalVaporFlux             = ff_param['nominalVaporFlux']
    outputsUpdate                = ff_param['outputsUpdate']
    bmapOutputUpdate             = ff_param['bmapOutputUpdate']

    # dimension
    nx = MNH_domain['nx']
    ny = MNH_domain['ny']
   
    #load attribute
    for key in ['Lx','Ly','Lz','NEx','NEy','NEz','SWx','SWy','SWz','t0','Lt']:
        attribute.domain[key] = MNH_domain[key]
    attribute.parameters['date'] =  MNH_domain['date']
    date_mnh = datetime.datetime.strptime(MNH_domain['date'], '%Y-%m-%d_%H:%M:%S')
    #time_mnh = attribute.domain['t0']
    attribute.parameters['refYear']  = date_mnh.year
    attribute.parameters['refMonth'] = date_mnh.month
    attribute.parameters['refDay']   = date_mnh.day
    
    #'projection','duration' are not defined here

    #add domain length to attribute and get atmo reso
    Lx = 1.*(attribute.domain['NEx']-attribute.domain['SWx'])
    Ly = 1.*(attribute.domain['NEy']-attribute.domain['SWy'])
    attribute.domain['Lx'] = Lx
    attribute.domain['Ly'] = Ly
    atmo_dx = Lx/nx
    atmo_dy = Ly/ny
    y_center = .5 * (Ly - 2 * atmo_dy) #of the domain without the grid pt on the side

    #init param
    windU_in = 10 # cte wind speed
    windV_in = 0
    mydatanc = './Inputs/mydata_'+expr_name+'.nc'
    mybmapnc = './Inputs/mybmap_'+expr_name+'.nc'

    #this need to be defined here for bmap matrix
    BMapsResolution = max(spatialIncrement/np.sqrt(2), minimalPropagativeFrontDepth)
    print 'BMapsResolution = ', BMapsResolution

    #---------------------------------------------
    #create data nc file for ForeFire Layer Models
    #---------------------------------------------
    print 'write data file: ', mydatanc.split('/')[-1]
    fnameout = mydatanc
    nc_dimension = {'DIMX': -999, 'DIMY': -999, 'DIMZ': -999, 'DIMT': -999, 'domdim': -999}


    var4dim = {'dim1': 'DIMX',   'dim2': 'DIMY',  'dim3': 'DIMZ', 'dim4': 'DIMT', 'type': 'float' }
    nc_variable  = {'vaporFlux' : var4dim, \
                    'heatFlux'  : var4dim, \
                    'FRP'       : var4dim, \
                    'fuel'      : var4dim, \
                    'windU'     : var4dim, \
                    'windV'     : var4dim, \
                    'altitude'  : var4dim, \
                    'parameters': {'dim1': 'domdim', 'type': 'c'},  
                    'domain'    : {'dim1': 'domdim', 'type': 'c'}  }

    nc_attribute  = {'vaporFlux': attribute.vaporFlux,
                    'heatFlux'  : attribute.heatFlux,\
                    'FRP'       : attribute.FRP,\
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
        elif key == 'heatFlux': 
            nco.variables[key][:] = np.zeros(dimensions_value) + 0
        elif key == 'FRP': 
            nco.variables[key][:] = np.zeros(dimensions_value) + 2
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
                values = np.zeros(dimensions_value[::-1]) -9999
                grid_y2d, grid_x2d = np.meshgrid(grid_y, grid_x)
                for i_line in range(nbre_fire_line):
                    x_line =  x_location_fireLine[i_line]
                    y_line =  length_fireLine[i_line]
                    idx_line_start = np.where( (grid_x2d>=(x_line  - depth_fireLine[i_line])) & (grid_x2d<= x_line             ) &\
                                               (grid_y2d>=(y_center-.5*y_line)              ) & (grid_y2d<=(y_center+.5*y_line))  )
                    try:
                        values[idx_line_start] = time_start_fireLine[i_line]
                        #plt.imshow(values.T,origin='lower',interpolation='nearest',\
                        #       extent=(grid_x.min(),grid_x.max()+BMapsResolution,grid_y.min(),grid_y.max()+BMapsResolution))
                        #plt.show()
                        #pdb.set_trace()
                    except IndexError:
                        pdb.set_trace()
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
    if flag_run_ff:
        ff = forefire.PLibForeFire()

        
        #set param
        ff.setString('ForeFireDataDirectory','Inputs')
        ff.setString('fireOutputDirectory','ForeFire/Outputs')
        ff.setInt('outputsUpdate',outputsUpdate)
        
        ff.setString('NetCDFfile',    mydatanc.split('/')[-1])
        ff.setString('fluxNetCDFfile',mydatanc.split('/')[-1])
        ff.setString('fuelsTableFile','fuels.ff')
        ff.setString('BMapFiles', mybmapnc.split('/')[-1])
        
        ff.setDouble("spatialIncrement",spatialIncrement)
        ff.setDouble("perimeterResolution",perimeterResolution)
        ff.setDouble("minimalPropagativeFrontDepth",minimalPropagativeFrontDepth)
        
        ff.setDouble("nominalHeatFlux",nominalHeatFlux)
        ff.setDouble("nominalVaporFlux",nominalVaporFlux)
        ff.setDouble("burningDuration",burningDuration)

        ff.setDouble("bmapOutputUpdate",bmapOutputUpdate)
        ff.setInt("defaultHeatType",0)
        ff.setInt("defaultFRPType",2)
        ff.setInt("defaultVaporType",3)

        #set domain
        ff.setInt("atmoNX",nx)
        ff.setInt("atmoNY",ny)
        ff.execute("FireDomain[sw=(%f,%f,0.);ne=(%f,%f,0.);t=%f]"%(attribute.domain['SWx'],attribute.domain['SWy'],\
                                                                    attribute.domain['NEx'],attribute.domain['NEy'],  \
                                                                    attribute.domain['t0']))

        extentLocal= getDomainExtent(ff.execute("print[]").split("\n")[0]);


        #set propagation model
        ff.addLayer("propagation","TroisPourcent","propagationModel")
        ff.addLayer("flux","heatFluxBasic","defaultHeatType")
        #ff.addLayer("flux","FRP","defaultFRPType")
        ff.addLayer("flux","vaporFluxBasic","defaultVaporType")

        
        fuelmap = ff.getDoubleArray("fuel").astype("int32")
        ff.addIndexLayer("table","fuel", extentLocal[0], extentLocal[2],0,  extentLocal[1]-extentLocal[0], extentLocal[3]-extentLocal[2], 0, fuelmap)

        print "resolution of bmap is ", ff.getString("bmapResolution")

        #set fire line
        if flag_extra_contour:
            print 'add an extra coutour'
            ff.execute("\tFireFront[t={%f}]".format(attribute.domain['t0']))
            ff.execute("\t\tFireNode[loc=(090,100,0.);vel=(0.,0.,0.);t=%f]".format(attribute.domain['t0']))
            ff.execute("\t\tFireNode[loc=(090,300,0.);vel=(0.,0.,0.);t=%f]".format(attribute.domain['t0']))
            ff.execute("\t\tFireNode[loc=(100,300,0.);vel=(0.,0.,0.);t=%f]".format(attribute.domain['t0']))
            ff.execute("\t\tFireNode[loc=(100,100,0.);vel=(0.,0.,0.);t=%f]".format(attribute.domain['t0']))

        #---------------------------------------------
        # run ForeFire simulation
        #---------------------------------------------
        pathes = []
        step = 10
        N_step = 20
        flux_out_ff_history = []
        for i in np.arange(1,N_step):
            
            ff_time = i*step

            print "goTo[t=%f]"%(i*step),
            ff.execute("goTo[t=%f]"%(i*step))
            
            #FRP = ff.getDoubleArray("FRP")[:,:,0]
            #plt.imshow(FRP.T,origin='lower',interpolation='nearest'); plt.show()

            flux2d = ff.getDoubleArray("heatFluxBasic")[:,:,0]
            flux_out_ff =  (flux2d * atmo_dx * atmo_dy).sum() * 1.e-6
            
            idx_fire_line = np.where(np.array(time_start_fireLine) <= ff_time)
            flux_expected = 0
            for idx_ in idx_fire_line[0]:
                flux_expected += depth_fireLine[idx_]*length_fireLine[idx_]*nominalHeatFlux * 1.e-6 
            print  '| flux out of ff = ', flux_out_ff, '| expected = ', flux_expected, '| ratio = ', flux_out_ff/flux_expected
            
            pathes += printToPathe( ff.execute("print[]"))
            flux_out_ff_history.append(flux_out_ff/flux_expected)

        ax = plt.subplot(111)
        ax.plot(step*np.arange(1,N_step),np.array(flux_out_ff_history))
        ax.set_ylabel('ratio:  heat_ff  /  heat_expected')
        ax.set_xlabel('time (s)')
        plt.show()
        sys.exit()
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


