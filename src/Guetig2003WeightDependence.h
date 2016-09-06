//
// Created by Simon Vogt on 05.09.16.
//

#ifndef AURYN_GUETIG2003WEIGHTDEPENDENCE_H
#define AURYN_GUETIG2003WEIGHTDEPENDENCE_H

#include "STDPWeightDependence.h"

namespace auryn
{
	/** \brief Scales synaptic weight updates by the exponential rule of Morrison2007 or even the Gütig2003 rule.
	 *
	 */
	class Guetig2003WeightDependence : public STDPWeightDependence
	{
		AurynFloat scale_LTP;
		AurynFloat scale_LTD;

		AurynFloat mu;

	public:
		Guetig2003WeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD, AurynFloat mu);
		virtual AurynDouble applyLTPscaling(AurynWeight *pDouble);
		virtual AurynDouble applyLTDscaling(AurynWeight *pDouble);
	};

}


#endif //AURYN_GUETIG2003WEIGHTDEPENDENCE_H
