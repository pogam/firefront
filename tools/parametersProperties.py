heatFlux = {'type':'flux',                   \
            'indices': [0],              \
            'model0name': 'heatFluxBasic',   \
            }
vaporFlux= {'type':'flux',                   \
            'indices': [3],              \
            'model3name': 'vaporFluxBasic',\
            }
fuel     = {'type' : "fuel"} 
windU    = {'type' : "data"} 
windV    = {'type' : "data"} 
altitude = {'type' : "data"} 

domain   = {'type': 'domain',
            'NEx'  : None ,
            'NEy'  : None, 
            'NEz'  : None, 
            'SWx' :  None,
            'SWy' :  None, 
            'SWz' :  None,
            'Lx' :  None,
            'Ly' :  None, 
            'Lz' :  None,
            't0'  : None, 
            'Lt'  : None 
            }

parameters= {'type' : "parameters" ,
            'date' : None  ,               #"2000-01-01_00:00:00",
            'refYear' : None,              #2000 ,
            'refMonth' : None,             #1,
            'refDay' : None,               #1,
            'duration' : None,             #360000 ,
	    'projection' : None,           #"OPENMAP" ,
            'projectionproperties' : None, #"41.551998,8.828396,41.551998,8.828396" ,
            }
