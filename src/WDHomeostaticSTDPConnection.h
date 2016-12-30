/*
 * Made by Simon
*/

#ifndef WDHOMEOSTATICSTDPCONNECTION_H_
#define WDHOMEOSTATICSTDPCONNECTION_H_

#include "auryn/auryn_definitions.h"
#include "auryn/AurynVector.h"
#include "auryn/DuplexConnection.h"
#include "auryn/Trace.h"
#include "auryn/LinearTrace.h"
#include "auryn/SpikeDelay.h"


#include "AdditiveWeightDependence.h"
#include "LinearAttractorWeightDependence.h"
#include "Guetig2003WeightDependence.h"
#include "Morrison2007WeightDependence.h"


namespace auryn {


/*! \brief STDP All-to-All Connection with pluggable weight dependence rules.
 *
 * This class implements standard STDP with a double exponential window and optinal
 * offset terms. Window amplitudes, time constants and weight dependence are freely configurable.
 */
class WDHomeostaticSTDPConnection : public DuplexConnection
{
public:
	class WeightDependentUpdateScaling;

private:
	//AurynFloat learningrate; /** Used to update the scaling factors in the weight dependence class */
	//AurynFloat A_minus; /*!< Amplitude of post-pre part of the STDP window */
	//AurynFloat A_plus; /*!< Amplitude of pre-post part of the STDP window */

protected:

	AurynFloat tau_pre;
	AurynFloat tau_post;
	Trace * tr_pre;
	Trace * tr_post;

	void propagate_forward();
	void propagate_backward();

	AurynWeight get_postspiketrace(NeuronID postspiker);
	AurynWeight get_prespiketrace(NeuronID prespiker);

	//const std::function<AurynDouble(AurynWeight)> &the_weight_dependence;
	const WeightDependentUpdateScaling* the_weight_dependence;


//	STDPWeightDependence* weight_dependence;
//	std::vector<HomeostasisRule*> homeostasisvector; ///< list of homeostasis rues to be applied to the connection in vector's order.

public:
	bool stdp_active;

	/** Creates a new connection with a default additive STDP rule. */
	WDHomeostaticSTDPConnection(SpikingGroup *source, NeuronGroup *destination);

	/** Creates a new connection with a given weight dependence rule. */
	WDHomeostaticSTDPConnection(SpikingGroup *source, NeuronGroup *destination, AurynWeight initialweight, AurynWeight maxweight,
								WeightDependentUpdateScaling *theWeightDependence, AurynFloat tau_pre=20e-3, AurynFloat tau_post=20e-3,
								AurynFloat sparseness=1.0, TransmitterType transmitter=GLUT);

	virtual ~WDHomeostaticSTDPConnection();
	virtual void finalize() override;
	void free();

	virtual void propagate() override;
	virtual void evolve() override;

	virtual string get_name() override;

	//virtual void set_max_weight(AurynWeight maximum_weight) override;

	void setTau_pre(AurynFloat tau_pre);
	void setTau_post(AurynFloat tau_post);

	//void evolve_homeostasis(HomeostasisRule* homeostasisobject);


		class WeightDependentUpdateScaling
		{
		private:

			// info for debug output and stuff:
			enum FactoryMethodUsed {Additive,Linear,Guetig2003,Morrison2007,Gilson2011};
			FactoryMethodUsed factory_used;

			WeightDependentUpdateScaling() {}	// never let anyone else call the default constructor. Use factory methods instead!
			//WeightDependentUpdateScaling(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize);

		protected:
			/* Constructors: */
			//WeightDependentUpdateScaling(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize);
			/** Protected Constructor: Objects of this type should only be produced via the factory methods provided below. */
			WeightDependentUpdateScaling(const AurynDouble (WeightDependentUpdateScaling::*pFunction_Causal)(const AurynWeight *)const,
										 const AurynDouble (WeightDependentUpdateScaling::*pFunction_Anticausal)(const AurynWeight *)const,
										 AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize);
			void compute_fudge(AurynWeight maxweight=1.0f);

			// caller-given parameters (used by every rule):
			AurynFloat A_plus;			// scale for causal spike order: pre-before-post
			AurynFloat A_minus;			// scale for anticausal spike order: post-before-pre
			AurynFloat update_stepsize; // learning rate

			// caller-given parameters (rule-specific):
			AurynFloat mu_1;
			AurynFloat mu_2;
			//AurynFloat attractorStrengthIndicator;	// don't really need to store these as I only use the fudge below.
			//AurynWeight attractorLocationIndicator;

			// fudge:
			AurynFloat scaleconstant_Causal;
			AurynFloat scaleconstant_Anticausal;
			AurynFloat slope_Causal;
			AurynFloat slope_Anticausal;
			AurynFloat offset_Causal;
			AurynFloat offset_Anticausal;


			/* Function pointers assigned by factory methods to avoid branching during execution of the simulation: */
			const AurynDouble (WeightDependentUpdateScaling::*scaleCausal)(const AurynWeight*)const; ///* Function pointer: takes AurynWeight and returns AurynDouble
			const AurynDouble (WeightDependentUpdateScaling::*scaleAnticausal)(const AurynWeight*)const; ///* Function pointer: takes AurynWeight and returns AurynDouble

			/* Scale pre-before-post updates: */
			const AurynDouble scaleCausal_Constant(    const AurynWeight* pWeight)const;
			const AurynDouble scaleCausal_Linear(      const AurynWeight* pWeight)const;
			const AurynDouble scaleCausal_Guetig2003(  const AurynWeight* pWeight)const;
			const AurynDouble scaleCausal_Morrison2007(const AurynWeight* pWeight)const;
			const AurynDouble scaleCausal_Gilson2011(  const AurynWeight* pWeight)const;

			/* Scale post-before-pre updates: */
			const AurynDouble scaleAnticausal_Constant(    const AurynWeight* pWeight)const;
			const AurynDouble scaleAnticausal_Linear(      const AurynWeight* pWeight)const;
			const AurynDouble scaleAnticausal_Guetig2003(  const AurynWeight* pWeight)const;
			const AurynDouble scaleAnticausal_Morrison2007(const AurynWeight* pWeight)const;
			const AurynDouble scaleAnticausal_Gilson2011(  const AurynWeight* pWeight)const;

		public:
			/* Factory methods to produce the weight updates you want: */
			static WeightDependentUpdateScaling* makeAdditiveUpdates(     AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize);
			static WeightDependentUpdateScaling* makeLinearUpdates(       AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat theAttractorStrengthIndicator, AurynWeight theAttractorLocationIndicator, AurynFloat theMeanSlope=0);
			static WeightDependentUpdateScaling* makeGuetig2003Updates(   AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat mu=0.02);
			static WeightDependentUpdateScaling* makeMorrison2007Updates( AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat mu_LTP=1, AurynFloat mu_LTD=0.2);
			static WeightDependentUpdateScaling* makeGilson2011Updates(   AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat mu_LTP, AurynFloat mu_LTD);

			const AurynDouble dothescale_Causal(const AurynWeight* pWeight)const;
			const AurynDouble dothescale_Anticausal(const AurynWeight* pWeight)const;

			void setMaxWeight(AurynWeight maxweight); ///* To set the max weight before giving this object to the Connection constractor, if needed.
			const string getRuleName()const;
		};

		class WeightDependentHomeostasisRule
		{

		};


	};

}

#endif /* WDHOMEOSTATICSTDPCONNECTION_H_ */
