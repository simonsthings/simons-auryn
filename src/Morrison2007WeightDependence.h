//
// Created by Simon Vogt on 05.09.16.
//

#ifndef AURYN_MORRISON2007WEIGHTDEPENDENCE_H
#define AURYN_MORRISON2007WEIGHTDEPENDENCE_H

#include "STDPWeightDependence.h"

namespace auryn
{
	/** \brief Scales synaptic weight updates by the exponential rule of Morrison2007 or even the Gütig2003 rule.
	 *
	 */
	class Morrison2007WeightDependence : public STDPWeightDependence
	{
		AurynFloat scale_LTP;
		AurynFloat scale_LTD;

		AurynFloat mu_LTP;
		AurynFloat mu_LTD;

	public:
		Morrison2007WeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD, AurynFloat mu_LTP, AurynFloat mu_LTD);
		virtual AurynDouble applyLTPscaling(AurynWeight *pDouble);
		virtual AurynDouble applyLTDscaling(AurynWeight *pDouble);
	};
}


#endif //AURYN_MORRISON2007WEIGHTDEPENDENCE_H
