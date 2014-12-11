/*

Copyright (C) 2012 ForeFire Team, SPE, Universit� de Corse.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 US

*/

#include "LavaPropagationModel.h"

namespace libforefire {

/* name of the model */
const string LavaPropagationModel::name = "LavaPropagationModel";

/* instantiation */
PropagationModel* getLavaPropagationModel(const int & mindex, DataBroker* db) {
	return new LavaPropagationModel(mindex, db);
}

/* registration */
int LavaPropagationModel::isInitialized =
        FireDomain::registerPropagationModelInstantiator(name, getLavaPropagationModel );

/* constructor */
LavaPropagationModel::LavaPropagationModel(const int & mindex, DataBroker* db)
: PropagationModel(mindex, db) {

	/* defining the properties needed for the model */
	effectiveSlope = registerProperty("slope");
//	heatFlux = registerProperty("heatFlux");
	viscosity = registerProperty("fuel.viscosity");
	flowSpeed = registerProperty("fuel.flowSpeed");
	rugosity = registerProperty("fuel.rugosity");
	curvature = registerProperty("frontCurvature");

	fastInSection = registerProperty("frontFastestInSection");

	/* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];

	/* registering the model in the data broker */
	dataBroker->registerPropagationModel(this);

	/* Definition of the coefficients */
	slowFactor = 1;
	lavaviscosity = 300;
	lavadensity = 3000;
	lavayield = 10000;
	if ( params->isValued("LavaModel.slowFactor") )
		slowFactor = params->getDouble("LavaModel.slowFactor");
/*	if ( params->isValued("lava.viscosity") )
		lavaviscosity = params->getDouble("lava.viscosity");
	if ( params->isValued("lava.density") )
		lavadensity = params->getDouble("lava.density");
	if ( params->isValued("lava.yield") )
		lavayield = params->getDouble("lava.yield");*/
}

/* destructor (shoudn't be modified) */
LavaPropagationModel::~LavaPropagationModel() {
}

/* accessor to the name of the model */
string LavaPropagationModel::getName(){
	return name;
}

/* *********************************************** */
/* Model for the propagation velocity of the front */
/* *********************************************** */

double LavaPropagationModel::getSpeed(double* valueOf){

//	if (valueOf[effectiveSlope] < 0){
//		return valueOf[viscosity]-valueOf[effectiveSlope]*valueOf[flowSpeed];
//	}
//	return valueOf[viscosity];
	double alpha = -atan(valueOf[effectiveSlope]);
//	double talpha = -valueOf[effectiveSlope];
//	double salpha = sin(alpha);
	double calpha= 1.57-alpha;
	double ctalpha= tan(calpha);
//	double csalpha= sin(calpha);
//	double alphadeg = alpha *180/3.14;
	double lavathick = (lavayield/(lavadensity*9.80*ctalpha));
	double speedLocale = (lavadensity*9.8*lavathick*lavathick*sin(calpha))/(3*lavaviscosity);
//	cout << "lavayield" << lavayield << "lavadensity" << lavadensity << endl;
//	cout <<"alpha "<< alpha << " alphadeg " << alphadeg << " lavathick " << lavathick << " speedLocale "<< speedLocale << endl;
	//cout<<valueOf[fastInSection]<<endl;


	if (valueOf[effectiveSlope] < 0){
		double speedL = speedLocale*valueOf[rugosity]+valueOf[flowSpeed];
		if (valueOf[fastInSection] > 0.99) return speedL;
	}
	return 0.01;

}

}

