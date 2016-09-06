//
// Created by Simon Vogt on 05.09.16.
//

#ifndef AURYN_LINEARWEIGHTDEPENDENCE_H
#define AURYN_LINEARWEIGHTDEPENDENCE_H

#include "STDPWeightDependence.h"

namespace auryn
{
	/** \brief Scales the STDP update by a linear function that depends on the current synaptic weight(s).
	 *
	 * TODO: reshape to allow vectorization?
	 */
	class LinearWeightDependence : public STDPWeightDependence
	{
		AurynFloat attractorStrengthIndicator;
		AurynWeight attractorLocationIndicator;

		AurynFloat fudge_LTP;
		AurynFloat fudge_LTD;

	public:
		LinearWeightDependence(AurynWeight w_max, AurynFloat theAttractorStrengthIndicator, AurynWeight theAttractorLocationIndicator);
		virtual AurynDouble applyLTPscaling(AurynWeight *pDouble);
		virtual AurynDouble applyLTDscaling(AurynWeight *pDouble);
	};

}


#endif //AURYN_LINEARWEIGHTDEPENDENCE_H
