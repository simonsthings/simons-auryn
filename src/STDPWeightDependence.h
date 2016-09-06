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
	private:

	protected:
		AurynWeight w_max;
		std::string rule_name;

	public:
		STDPWeightDependence();
		STDPWeightDependence(AurynWeight w_max);
		virtual AurynDouble applyLTPscaling(AurynWeight *pDouble)=0;
		virtual AurynDouble applyLTDscaling(AurynWeight *pDouble)=0;
		AurynWeight getW_max() const;
	};

}




#endif //AURYN_STDPWEIGHTDEPENDENCE_H
