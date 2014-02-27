/*

Copyright (C) 2012 ForeFire Team, SPE, Université de Corse.

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

#ifndef LAVAHEATFLUXMODEL_H_
#define LAVAHEATFLUXMODEL_H_

#include "FluxModel.h"
#include "FireDomain.h"
#include "Futils.h"

namespace libforefire {

class LavaHeatFluxModel: public FluxModel {

	/*! name the model */
	static const string name;

	/*! boolean for initialization */
	static int isInitialized;

	/*! properties needed by the model */
	size_t windU;
	size_t windV;

	/*! coefficients needed by the model */
	double eruptionTime;
	double crustTemperature;
	double lavaTemperature;
	vector<double> refHours;
	vector<double> crustFractions;
	vector<double> windValues;
	vector<double> A;
	vector<double> B;

	/*! local variables */

	/*! result of the model */
	double getValue(double*, const double&
			, const double&, const double&);

public:

	LavaHeatFluxModel(const int& = 0, DataBroker* = 0);
	virtual ~LavaHeatFluxModel();

	string getName();
};

FluxModel* getLavaHeatFluxModel(const int& = 0, DataBroker* = 0);

} /* namespace libforefire */
#endif /* LAVAHEATFLUXMODEL_H_ */
