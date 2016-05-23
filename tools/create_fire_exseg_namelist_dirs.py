import os
import sys
import f90nml

##################################################
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)



###########################################
if __name__ == '__main__':
###########################################
   
    #create ForeFire dir
    ensure_dir('./ForeFire/')
   
    #read MNH info
    #-------------
    MNH_namelist = f90nml.read('../02_mnh/EXSEG1.nam')
    expr_name = MNH_namelist['NAM_CONF']['CEXP']

    #read extra namelist for ff coupling
    #-------------
    MNH_Fire_namelist_extra = f90nml.read('./ff-mnh.nam')

    #update and add keys
    #-------------
    keys = list(set(MNH_namelist.keys()+MNH_Fire_namelist_extra.keys()))
    for key in keys:
        if key in MNH_Fire_namelist_extra.keys() :
            MNH_namelist[key] = MNH_Fire_namelist_extra[key]
    
    #create exseg file for mnh with ff
    #-------------
    MNH_namelist.write('./EXSEG1.nam',force=True)

    #read FF info
    #-------------
    nmls = f90nml.read('../03_ff/Inputs/ff-param.nam')
    ff_param = nmls['FF_PARAM']


    #create FF environment
    #-------------
    ForeFireDataDirectory = 'ForeFire'
    key_needed_in_Paramsff = ['ForeFireDataDirectory','fireOutputDirectory','atmoOutputDirectories',\
                              'outputsUpdate','NetCDFfile','fuelsTableFile','fluxNetCDFfile','BMapFiles',\
                              'bmapOutputUpdate','propagationModel','spatialIncrement','perimeterResolution',\
                              'minimalPropagativeFrontDepth','burningDuration','nominalHeatFlux',\
                              'nominalVaporFlux','watchedProc','outputDirectories']
    
    #write Params.ff
    def print_param_in_paramsFF_file(f,key,value):
        f.write('setParameters[{:s}={:s}]\n'.format(key,value) )

    f=open('./ForeFire/Params.ff','w')
    for key in key_needed_in_Paramsff:
        if key == 'ForeFireDataDirectory':
            print_param_in_paramsFF_file(f,key,ForeFireDataDirectory)
        elif key == 'NetCDFfile':
            print_param_in_paramsFF_file(f,key,'mydata_'+expr_name+'.nc')
            if os.path.isfile(ForeFireDataDirectory+'/mydata_'+expr_name+'.nc'):
                os.remove(ForeFireDataDirectory+'/mydata_'+expr_name+'.nc')
            os.symlink('../../03_ff/Inputs/mydata_'+expr_name+'.nc',ForeFireDataDirectory+'/mydata_'+expr_name+'.nc')
        elif key == 'fluxNetCDFfile':
            print_param_in_paramsFF_file(f,key,'mydata_'+expr_name+'.nc')
        elif key == 'BMapFiles':
            print_param_in_paramsFF_file(f,key,'mybmap_'+expr_name+'.nc')
            if os.path.isfile(ForeFireDataDirectory+'/mybmap_'+expr_name+'.nc'):
                os.remove(ForeFireDataDirectory+'/mybmap_'+expr_name+'.nc')
            os.symlink('../../03_ff/Inputs/mybmap_'+expr_name+'.nc',ForeFireDataDirectory+'/mybmap_'+expr_name+'.nc')
        else:
            print_param_in_paramsFF_file(f,key,str(ff_param[key]))
            if key == 'fuelsTableFile':
                if os.path.isfile(ForeFireDataDirectory+'/'+ff_param[key]):
                    os.remove(ForeFireDataDirectory+'/'+ff_param[key])
                os.symlink('../../03_ff/Inputs/'+ff_param[key],ForeFireDataDirectory+'/'+ff_param[key])
            elif key in ['fireOutputDirectory','atmoOutputDirectories','outputDirectories']:
                ensure_dir(ff_param[key]+'/')

    f.close()
    
    #write Init.ff
    f=open('./ForeFire/Init.ff','w')
    f.close()
    
    
