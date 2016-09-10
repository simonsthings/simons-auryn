//
// Created by Simon Vogt on 05.09.16.
//

#include "Morrison2007WeightDependence.h"

using namespace auryn;


/* Morrison 2007 weight dependence: */
Morrison2007WeightDependence::Morrison2007WeightDependence(AurynFloat mu_LTP, AurynFloat mu_LTD)
		: mu_LTP(mu_LTP), mu_LTD(mu_LTD)
{}
Morrison2007WeightDependence::Morrison2007WeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD, AurynFloat mu_LTP, AurynFloat mu_LTD)
		: STDPWeightDependence(w_max,scale_LTP,scale_LTD), mu_LTP(mu_LTP), mu_LTD(mu_LTD)
{}
AurynDouble Morrison2007WeightDependence::scalePreBeforePost(AurynWeight *weight)
{
	return std::pow(*weight,mu_LTP) * scaleconstant_PreBeforePost ;
}
AurynDouble Morrison2007WeightDependence::scalePreAfterPost(AurynWeight *weight)
{
	return std::pow(*weight,mu_LTD) * scaleconstant_PreAfterPost ;
}





