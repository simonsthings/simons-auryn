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
	class LinearAttractorWeightDependence : public STDPWeightDependence
	{
	private:
		AurynFloat attractorStrengthIndicator;
		AurynWeight attractorLocationIndicator;

		AurynFloat fudge_LTP;
		AurynFloat fudge_LTD;

		void compute_fudge();

	public:
		LinearAttractorWeightDependence(AurynFloat theAttractorStrengthIndicator=0.5, AurynWeight theAttractorLocationIndicator=0.5);
		LinearAttractorWeightDependence(AurynWeight w_max, AurynFloat scale_LTP, AurynFloat scale_LTD, AurynFloat theAttractorStrengthIndicator, AurynWeight theAttractorLocationIndicator);
		virtual AurynDouble scalePreBeforePost(AurynWeight* weight);
		virtual AurynDouble scalePreAfterPost(AurynWeight* weight);
	};

}


#endif //AURYN_LINEARWEIGHTDEPENDENCE_H
