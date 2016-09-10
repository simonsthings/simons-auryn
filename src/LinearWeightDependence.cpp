//
// Created by Simon Vogt on 05.09.16.
//

#include "LinearWeightDependence.h"

using namespace auryn;


/* Linear weight dependence with given slopes: */
LinearWeightDependence::LinearWeightDependence(AurynFloat theAttractorStrengthIndicator, AurynWeight theAttractorLocationIndicator)
		: attractorStrengthIndicator(theAttractorStrengthIndicator), attractorLocationIndicator(theAttractorLocationIndicator)
{
	compute_fudge();
}
LinearWeightDependence::LinearWeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD, AurynFloat theAttractorStrengthIndicator, AurynWeight theAttractorLocationIndicator)
		: STDPWeightDependence(w_max,scale_LTP,scale_LTD), attractorStrengthIndicator(theAttractorStrengthIndicator), attractorLocationIndicator(theAttractorLocationIndicator)
{
	compute_fudge();
}
void LinearWeightDependence::compute_fudge()
{
	fudge_LTP = 0.5f*w_max-(w_max-attractorLocationIndicator)*attractorStrengthIndicator;
	fudge_LTD = 0.5f*w_max-(attractorLocationIndicator)*attractorStrengthIndicator;
}
AurynDouble LinearWeightDependence::scalePreBeforePost(AurynWeight *weight)
{
	return scaleconstant_PreBeforePost * ( (w_max - *weight) * attractorStrengthIndicator + fudge_LTP );
}
AurynDouble LinearWeightDependence::scalePreAfterPost(AurynWeight *weight)
{
	return scaleconstant_PreAfterPost  * ( (*weight) * attractorStrengthIndicator + fudge_LTD );
}



