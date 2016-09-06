//
// Created by Simon Vogt on 05.09.16.
//

#ifndef AURYN_ADDITIVEWEIGHTDEPENDENCE_H
#define AURYN_ADDITIVEWEIGHTDEPENDENCE_H

#include "STDPWeightDependence.h"

namespace auryn
{
	/** \brief Scales synaptic weight updates by the exponential rule of Morrison2007 or even the Gütig2003 rule.
	 *
	 */
	class AdditiveWeightDependence : public STDPWeightDependence
	{
		AurynFloat scale_LTP;
		AurynFloat scale_LTD;

	public:
		AdditiveWeightDependence(AurynWeight w_max, AurynFloat scale_LTP=1, AurynFloat scale_LTD=1);
		virtual AurynDouble applyLTPscaling(AurynWeight *pDouble);
		virtual AurynDouble applyLTDscaling(AurynWeight *pDouble);
	};
}


#endif //AURYN_ADDITIVEWEIGHTDEPENDENCE_H
