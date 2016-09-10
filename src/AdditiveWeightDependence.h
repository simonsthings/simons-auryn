//
// Created by Simon Vogt on 05.09.16.
//

#ifndef AURYN_ADDITIVEWEIGHTDEPENDENCE_H
#define AURYN_ADDITIVEWEIGHTDEPENDENCE_H

#include "STDPWeightDependence.h"

namespace auryn
{
	/** \brief This rule produces additive STDP, so there is actually no dependence on the given weight.
	 *
	 * Only the given scale constants for each side of the STDP rule are returned by calls to scalePreBeforePost() and scalePreAfterPost().
	 */
	class AdditiveWeightDependence : public STDPWeightDependence
	{

	public:
		AdditiveWeightDependence(AurynWeight w_max=1, AurynFloat scale_LTP=1, AurynFloat scale_LTD=1);
		virtual AurynDouble scalePreBeforePost(AurynWeight *pDouble);
		virtual AurynDouble scalePreAfterPost(AurynWeight *pDouble);
	};
}


#endif //AURYN_ADDITIVEWEIGHTDEPENDENCE_H
