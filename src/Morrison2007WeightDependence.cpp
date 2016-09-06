//
// Created by Simon Vogt on 05.09.16.
//

#include "Morrison2007WeightDependence.h"

using namespace auryn;


/* Morrison 2007 weight dependence: */
Morrison2007WeightDependence::Morrison2007WeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD, AurynFloat mu_LTP, AurynFloat mu_LTD)
		: STDPWeightDependence(w_max), scale_LTP(scale_LTP), scale_LTD(scale_LTD), mu_LTP(mu_LTP), mu_LTD(mu_LTD)
{}
AurynDouble Morrison2007WeightDependence::applyLTPscaling(AurynWeight *weight)
{
	return std::pow(*weight,mu_LTP) * scale_LTP ;
}
AurynDouble Morrison2007WeightDependence::applyLTDscaling(AurynWeight *weight)
{
	return std::pow(*weight,mu_LTD) * scale_LTD ;
}




