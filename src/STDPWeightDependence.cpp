//
// Created by Simon Vogt on 05.09.16.
//

#include "STDPWeightDependence.h"

using namespace auryn;

STDPWeightDependence::STDPWeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD)
		: w_max(w_max), scaleconstant_PreBeforePost(scale_LTP), scaleconstant_PreAfterPost(scale_LTD)
{
}

STDPWeightDependence::~STDPWeightDependence()
{
}











