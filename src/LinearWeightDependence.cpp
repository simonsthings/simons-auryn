//
// Created by Simon Vogt on 05.09.16.
//

#include "LinearWeightDependence.h"

using namespace auryn;


/* Linear weight dependence with given slopes: */
LinearWeightDependence::LinearWeightDependence(AurynWeight w_max, AurynFloat theAttractorStrengthIndicator, AurynWeight theAttractorLocationIndicator)
		: STDPWeightDependence(w_max), attractorStrengthIndicator(theAttractorStrengthIndicator), attractorLocationIndicator(theAttractorLocationIndicator)
{
	fudge_LTP = attractorStrengthIndicator+0.5f*w_max-(w_max-attractorLocationIndicator)*attractorStrengthIndicator;
	fudge_LTD = attractorStrengthIndicator+0.5f*w_max-(attractorLocationIndicator)*attractorStrengthIndicator;
}
AurynDouble LinearWeightDependence::applyLTPscaling(AurynWeight *weight)
{
	return (w_max - *weight) * fudge_LTP ;
}
AurynDouble LinearWeightDependence::applyLTDscaling(AurynWeight *weight)
{
	return (*weight) * fudge_LTD ;
}

