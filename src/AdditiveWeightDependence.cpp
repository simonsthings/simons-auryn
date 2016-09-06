//
// Created by Simon Vogt on 05.09.16.
//

#include "AdditiveWeightDependence.h"

using namespace auryn;

AdditiveWeightDependence::AdditiveWeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD)
		: STDPWeightDependence(w_max), scale_LTP(scale_LTP), scale_LTD(scale_LTD)
{}
AurynDouble AdditiveWeightDependence::applyLTPscaling(AurynWeight *weight)
{
	return scale_LTP ;
}
AurynDouble AdditiveWeightDependence::applyLTDscaling(AurynWeight *weight)
{
	return scale_LTD ;
}
