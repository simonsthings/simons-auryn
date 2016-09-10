//
// Created by Simon Vogt on 05.09.16.
//

#include "Guetig2003WeightDependence.h"

using namespace auryn;


/* Gütig 2003 weight dependence: */
Guetig2003WeightDependence::Guetig2003WeightDependence(AurynFloat mu) : mu(mu)
{}
Guetig2003WeightDependence::Guetig2003WeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD, AurynFloat mu)
		: STDPWeightDependence(w_max,scale_LTP,scale_LTD), mu(mu)
{}
AurynDouble Guetig2003WeightDependence::scalePreBeforePost(AurynWeight *weight)
{
	return std::pow(1-*weight,mu) * scaleconstant_PreBeforePost ;
}
AurynDouble Guetig2003WeightDependence::scalePreAfterPost(AurynWeight *weight)
{
	return std::pow(*weight,mu) * scaleconstant_PreAfterPost ;
}



