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

#include "ForeFireV1HeatFluxModel.h"

namespace libforefire {

/* name of the model */
const string ForeFireV1HeatFluxModel::name = "ForeFireV1HeatFlux";

/* registration */
int ForeFireV1HeatFluxModel::isInitialized =
        FireDomain::registerFluxModelInstantiator(name, getForeFireV1HeatFluxModel );

/* instantiation */
FluxModel* getForeFireV1HeatFluxModel(const int& index, DataBroker* db) {
	return new ForeFireV1HeatFluxModel(index, db);
}

/* constructor */
ForeFireV1HeatFluxModel::ForeFireV1HeatFluxModel(
		const int & mindex, DataBroker* db) : FluxModel(mindex, db) {

	/* defining the properties needed for the model */
	sd = registerProperty("fuel.sd");
	md = registerProperty("fuel.Md");
	ml = registerProperty("fuel.Ml");
	sigmad = registerProperty("fuel.Sigmad");
	sigmal = registerProperty("fuel.Sigmal");
	lai = registerProperty("fuel.Blai");
	DeltaH = registerProperty("fuel.DeltaH");
	m_e = registerProperty("fuel.me");

	/* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];

	/* registering the model in the data broker */
	dataBroker->registerFluxModel(this);

	/* Definition of the coefficients */
	BurrowsCst = 208487;
	if ( params->isValued("FFfluxes.BurrowsCst") )
		BurrowsCst = params->getDouble("FFfluxes.BurrowsCst");
	liveTimeScaleRatio = 5.;
	if ( params->isValued("FFfluxes.liveTimeScaleRatio") )
		liveTimeScaleRatio = params->getDouble("FFfluxes.liveTimeScaleRatio");
	if ( !params->isValued("FFfluxes.LAICoeffs") )
		cout<<"ERROR: vector of parameters FFfluxes.LAICoeffs should be valued"<<endl;
	LAICoeffs = params->getDoubleArray("FFfluxes.LAICoeffs");
	if ( !params->isValued("FFfluxes.mCoeffs") )
		cout<<"ERROR: vector of parameters FFfluxes.mCoeffs should be valued"<<endl;
	mCoeffs = params->getDoubleArray("FFfluxes.mCoeffs");
	chi_r = 0.3;
	if ( params->isValued("radiativeFraction") )
		chi_r = params->getDouble("radiativeFraction");
	chi_b = 1.;
	if ( params->isValued("FFfluxes.chi_b") )
		chi_b = params->getDouble("FFfluxes.chi_b");
	 cout << "new V1 Heat Fl Model "<< endl;
}

/* destructor (shoudn't be modified) */
ForeFireV1HeatFluxModel::~ForeFireV1HeatFluxModel() {
	if ( properties != 0 ) delete properties;
}

/* accessor to the name of the model */
string ForeFireV1HeatFluxModel::getName(){
	return name;
}

double ForeFireV1HeatFluxModel::texp(double& t, double& tau){
	return t/tau/tau*exp(-t/tau);
}

double ForeFireV1HeatFluxModel::tauexp(double& t, double& tau){
	return (t+tau)/tau*exp(-t/tau);
}

/* ****************** */
/* Model for the flux */
/* ****************** */

double ForeFireV1HeatFluxModel::getValue(double* valueOf
		, const double& bt, const double& et, const double& at){

	/* see ForeFire manuals for flux models */

	double T0 = BurrowsCst/pow(valueOf[sd],1.236);
	double LAICorr = 1. + LAICoeffs[0]*(valueOf[lai]-2)*(valueOf[lai]-2)
			/(valueOf[lai]+LAICoeffs[1])/(valueOf[lai]+LAICoeffs[1]);
	double moistCorr = mCoeffs[0] + mCoeffs[1]*exp(mCoeffs[2]*valueOf[md])
		 + mCoeffs[3]*exp(-mCoeffs[4]*valueOf[md]);
	double taud = LAICorr*moistCorr*T0;
	double taul = liveTimeScaleRatio*taud;
	double deadPotHeat = (1.-chi_r)*valueOf[sigmad]*valueOf[DeltaH];
	double livePotHeat = chi_b*exp(-4.*valueOf[ml]/m_e)*valueOf[sigmal]*valueOf[DeltaH];

	/* Instantaneous flux */
	/* ------------------ */
	if ( bt == et ){
		if ( bt < at ) return 0;
		double dtb = bt - at;
		return deadPotHeat*texp(dtb,taud) + livePotHeat*texp(dtb,taul);
	}

	/* Averaged flux */
	/* ------------- */
	if ( et < at ) return 0;
	double dtb = bt - at;
	if ( dtb < 0 ) dtb = 0.;
	double dte = et - at;
	return ( deadPotHeat*(tauexp(dtb,taud)-tauexp(dte,taud))
		   + livePotHeat*(tauexp(dtb,taul)-tauexp(dte,taul)) )/(et-bt);
}

} /* namespace libforefire */
