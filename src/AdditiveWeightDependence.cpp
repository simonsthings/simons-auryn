//
// Created by Simon Vogt on 05.09.16.
//

#include "AdditiveWeightDependence.h"

using namespace auryn;

AdditiveWeightDependence::AdditiveWeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD)
		: STDPWeightDependence(w_max,scale_LTP,scale_LTD)
{}
AurynDouble AdditiveWeightDependence::scalePreBeforePost(AurynWeight *weight)
{
	// ignore the given weight:
	return scaleconstant_PreBeforePost;
}
AurynDouble AdditiveWeightDependence::scalePreAfterPost(AurynWeight *weight)
{
	// ignore the given weight:
	return scaleconstant_PreAfterPost;
}
