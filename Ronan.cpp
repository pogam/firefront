/*

Copyright (C) 2012 ForeFire Team, SPE, UniversitŽ de Corse.

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

#include "Ronan.h"

using namespace std;

namespace libforefire {

/* name of the model */
const string Ronan::name = "Ronan";

/* instantiation */
PropagationModel* getRonanModel(const int & mindex, DataBroker* db) {
	return new Ronan(mindex, db);
}

/* registration */
int Ronan::isInitialized =
		FireDomain::registerPropagationModelInstantiator(name, getRonanModel );

/* constructor */
Ronan::Ronan(const int & mindex, DataBroker* db)
: PropagationModel(mindex, db) {
	/* defining the properties needed for the model */

	/* allocating the vector for the values of these properties */
	if ( numProperties > 0 ) properties =  new double[numProperties];

	/* registering the model in the data broker */
	dataBroker->registerPropagationModel(this);

	/* Definition of the coefficients */
}

/* destructor (shoudn't be modified) */
Ronan::~Ronan() {
}

/* accessor to the name of the model */
string Ronan::getName(){
	return name;
}

/* *********************************************** */
/* Model for the propagation velovity of the front */
/* *********************************************** */

double Ronan::getSpeed(double* valueOf){

	return 0;

}

} /* namespace libforefire */
