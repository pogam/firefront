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

#include "LavaHCLFluxModel.h"

using namespace std;

namespace libforefire {

/* name of the model */
const string LavaHCLFluxModel::name = "LavaHCLFlux";

/* instantiation */
FluxModel* getLavaHCLFluxModel(const int& index, DataBroker* db) {
	return new LavaHCLFluxModel(index, db);
}

/* registration */
int LavaHCLFluxModel::isInitialized =
        FireDomain::registerFluxModelInstantiator(name, getLavaHCLFluxModel );

/* constructor */
LavaHCLFluxModel::LavaHCLFluxModel(
		const int & mindex, DataBroker* db)
	: FluxModel(mindex, db) {
	/* defining the properties needed for the model */
	/* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];
	/* registering the model in the data broker */
	dataBroker->registerFluxModel(this);

	/* Definition of the coefficients */
	arrivalTime = 5600.;
	if(params->isValued("vaporFlux.activeArea"))
//		cout<<params->getDouble("vaporFlux.activeArea")<<endl;
	if ( params->isValued("laze.arrivalTime") )
		arrivalTime = params->getDouble("laze.arrivalTime");
	if ( !params->isValued("laze.hours") )
		cout<<"ERROR: vector of parameters laze.hours should be valued"<<endl;
	refHours = params->getDoubleArray("laze.hours");
	if ( !params->isValued("laze.flows") )
		cout<<"ERROR: vector of parameters laze.flows should be valued"<<endl;
	refFlows = params->getDoubleArray("laze.flows");
	exchangeArea = 80000.;
	if ( params->isValued("laze.exchangeArea") )
		exchangeArea = params->getDouble("laze.exchangeArea");
	/* coefficients */
	// Dividing by the final area of exchange to convert to kg.m-2.s-1
    convert= 1./exchangeArea;
}

/* destructor (shoudn't be modified) */
LavaHCLFluxModel::~LavaHCLFluxModel() {
	if ( properties != 0 ) delete properties;
}

/* accessor to the name of the model */
string LavaHCLFluxModel::getName(){
	return name;
}

/* ****************** */
/* Model for the flux */
/* ****************** */

double LavaHCLFluxModel::getValue(double* valueOf
		, const double& bt, const double& et, const double& at){

	if ( bt - arrivalTime < 0 ) return 0.;


	/* getting the hours since eruption */
	double hoursSinceArrival = (bt-arrivalTime)/3600.;
	if ( hoursSinceArrival > refHours.back() ) return 0.;

	/* getting the index in the vector of fluxes */
	size_t hind = 0;
	size_t nhours = refHours.size();
	while ( hind+1 < nhours and refHours[hind+1] < hoursSinceArrival ) hind++;
	/* interpolation of the flux between the values */
	double beta = (hoursSinceArrival-refHours[hind])
			/(refHours[hind+1]-refHours[hind]);
	double flux = beta*refFlows[hind+1] + (1.-beta)*refFlows[hind];
	double HCL = flux * 0.00033 *0.5; 
//	cout << " fluxhcl" << flux << endl;
//	cout << "LavaLaze " << params->getDouble("LavaHCLFlux.activeArea") +1 << endl ;
	if((bt-at) < (35*3600)) return HCL /((params->getDouble("LavaHCLFlux.activeArea"))+1.);
	       return 0;
//	return convert*flux/vaporFlux.activeArea;

}

} /* namespace libforefire */
