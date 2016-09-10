//
// Created by Simon Vogt on 05.09.16.
//

#ifndef AURYN_GUETIG2003WEIGHTDEPENDENCE_H
#define AURYN_GUETIG2003WEIGHTDEPENDENCE_H

#include "STDPWeightDependence.h"

namespace auryn
{
	/** \brief Scales synaptic weight updates by the exponential rule of Gütig2003.
	 *
	 */
	class Guetig2003WeightDependence : public STDPWeightDependence
	{
		AurynFloat mu;

	public:
		Guetig2003WeightDependence(AurynFloat mu=0.02);
		Guetig2003WeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD, AurynFloat mu);
		virtual AurynDouble scalePreBeforePost(AurynWeight *pDouble);
		virtual AurynDouble scalePreAfterPost(AurynWeight *pDouble);
	};

}


#endif //AURYN_GUETIG2003WEIGHTDEPENDENCE_H
