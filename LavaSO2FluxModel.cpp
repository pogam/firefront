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

#include "LavaSO2FluxModel.h"

using namespace std;

namespace libforefire {

/* name of the model */
const string LavaSO2FluxModel::name = "LavaSO2Flux";

/* instantiation */
FluxModel* getLavaSO2FluxModel(const int& index, DataBroker* db) {
	return new LavaSO2FluxModel(index, db);
}

/* registration */
int LavaSO2FluxModel::isInitialized =
        FireDomain::registerFluxModelInstantiator(name, getLavaSO2FluxModel );

/* constructor */
LavaSO2FluxModel::LavaSO2FluxModel(
		const int & mindex, DataBroker* db)
	: FluxModel(mindex, db) {
	/* defining the properties needed for the model */
	/* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];
	/* registering the model in the data broker */
	dataBroker->registerFluxModel(this);

	/* Definition of the coefficients */
	eruptionTime = 0.;
	if ( params->isValued("lava.eruptionTime") )
		eruptionTime = params->getDouble("lava.eruptionTime");
	lavaArea = 5.5e6;
	if ( params->isValued("lava.area") )
		lavaArea = params->getDouble("lava.area");
	if ( !params->isValued("SO2.hours") )
		cout<<"ERROR: vector of parameters SO2.hours should be valued"<<endl;
	refHours = params->getDoubleArray("SO2.hours");
	if ( !params->isValued("SO2.flows") )
		cout<<"ERROR: vector of parameters SO2.flows should be valued"<<endl;
	refFlows = params->getDoubleArray("SO2.flows");
	emissionRatio = 0.;
	if ( params->isValued("SO2.craterRatio") )
		emissionRatio = 1. - params->getDouble("SO2.craterRatio");

	/* local variables */
	// Error estimate correction
	convert = 392./231.;
	// Dividing by the final area of the lava to convert to kg.m-2.s-1
	convert = convert/lavaArea;
	// converting from kg.m-2.s-1 to molecules.m-2.s-1
//	convert = convert/64.e-3;
//	convert = convert/64.e-3*6.022e23;
}

/* destructor (shoudn't be modified) */
LavaSO2FluxModel::~LavaSO2FluxModel() {
	if ( properties != 0 ) delete properties;
}

/* accessor to the name of the model */
string LavaSO2FluxModel::getName(){
	return name;
}

/* ****************** */
/* Model for the flux */
/* ****************** */

double LavaSO2FluxModel::getValue(double* valueOf
		, const double& bt, const double& et, const double& at){
	if ( bt - eruptionTime < 0 ) return 0.;

	/* getting the hours since eruption */
	double hoursSinceEruption = (bt-eruptionTime)/3600.;
	if ( hoursSinceEruption > refHours.back() ) return 0.;

	/* getting the index in the vector of fluxes */
	size_t hind = 0;
	size_t nhours = refHours.size();
	while ( hind+1 < nhours and refHours[hind+1] < hoursSinceEruption ) hind++;
	/* interpolation of the flux between the values */
//	double beta = (hoursSinceEruption-refHours[hind])
//			/(refHours[hind+1]-refHours[hind]);
	//double flux = beta*refFlows[hind+1] + (1.-beta)*refFlows[hind];
//	cout << "lavaso2" <<  lavaso2  << endl;
	return 0;

}

} /* namespace libforefire */
