//
// Created by Simon Vogt on 05.09.16.
//

#ifndef AURYN_STDPWEIGHTDEPENDENCE_H
#define AURYN_STDPWEIGHTDEPENDENCE_H

#include <auryn/auryn_definitions.h>

namespace auryn
{

	/*! \brief Generalized weight dependence base class for STDP rules.
	 *
	 * This class implements a number of weight dependencies that are used by the GeneralAlltoallSTDPConnection class.
	 * This allows us to mix-n-match STDP rules with different implementations of weight dependence.
	 */
	class STDPWeightDependence
	{
	public:
		/// Public fields (because why not):
		AurynFloat scaleconstant_PreBeforePost;
		AurynFloat scaleconstant_PreAfterPost;

		AurynWeight w_max;
		std::string rule_name;

		/// Constructors:
		STDPWeightDependence(AurynWeight w_max=1, AurynFloat scale_LTP=1, AurynFloat scale_LTD=1);
		virtual ~STDPWeightDependence();
		virtual AurynDouble scalePreBeforePost(AurynWeight *pDouble)=0;
		virtual AurynDouble scalePreAfterPost(AurynWeight *pDouble)=0;
	};

}

#endif //AURYN_STDPWEIGHTDEPENDENCE_H

