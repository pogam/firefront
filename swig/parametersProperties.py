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
            'NEx'  : 420.0,
            'NEy'  : 420.0, 
            'NEz'  :    0.0, 
            'SWx' :  -20.0,
            'SWy' :  -20.0, 
            'SWz' :   0.0,
            'Lx' :  None,
            'Ly' :  None, 
            'Lz' :  None,
            't0'  :   0.0, 
            'Lt'  : 1000.0, 
            }
parameters= {'type' : "parameters" ,
            'date' : "2000-01-01_00:00:00",
            'refYear' : 2000 ,
            'refDay' : 0,
            'duration' : 360000 ,
	    'projection' : "OPENMAP" ,
            'projectionproperties' : "41.551998,8.828396,41.551998,8.828396" ,
            }
