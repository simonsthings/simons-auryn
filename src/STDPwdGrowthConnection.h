/*
 * Made by Simon
*/

#ifndef STDPWDHOMCONNECTION_H_
#define STDPWDHOMCONNECTION_H_

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


namespace auryn
{
	/*! \brief STDP All-to-All Connection with pluggable weight dependence rules.
	 *
	 * This class implements standard STDP with a double exponential window and optinal
	 * offset terms. Window amplitudes, time constants and weight dependence are freely configurable.
	 */
	class STDPwdGrowthConnection : public DuplexConnection
	{
	public:
		/* Forward declarations; see below */
		class WeightDependentUpdatescalingRule;
		class WeightDependentGrowthRule;

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
		WeightDependentUpdatescalingRule* the_weight_dependence; // required in constructor.

	//	std::vector<HomeostasisRule*> homeostasisvector; ///< list of homeostasis rules to be applied to the connection in vector's order.
		WeightDependentGrowthRule* the_growth_rule = nullptr;  // init with NULL

		/** A one-bit-per-bool storage for remembering whether any incoming weights of each postsynaptic neuron have changed.
		 * This allows us to only perform the computeTrainedness() of the GrowthRule on that subset of postsynaptic neurons for speed gains.
		 * (Updates during STDP deactivated for now because possibly not needed anyway..) */
		std::vector<bool> postsyn_weights_have_changed;

	public:
		bool stdp_active;

		/** Creates a new connection with a default additive STDP rule. */
		STDPwdGrowthConnection(SpikingGroup *source, NeuronGroup *destination);

		/** Creates a new connection with a given weight dependence rule. */
		STDPwdGrowthConnection(SpikingGroup *source, NeuronGroup *destination, AurynWeight initialweight, AurynWeight maxweight,
									WeightDependentUpdatescalingRule *theWeightDependence, AurynFloat tau_pre=20e-3, AurynFloat tau_post=20e-3,
									AurynFloat sparseness=1.0, TransmitterType transmitter=GLUT);

		virtual ~STDPwdGrowthConnection();
		virtual void finalize() override;

		virtual void propagate() override;
		virtual void evolve() override;
		virtual string get_name() override;

		void setUpdateScalingRule(WeightDependentUpdatescalingRule* pUpdatescalingRule);
		void setGrowthRule(WeightDependentGrowthRule* pGrowthRule);

		//virtual void set_max_weight(AurynWeight maximum_weight) override;

		void setTau_pre(AurynFloat tau_pre);
		void setTau_post(AurynFloat tau_post);

		//void evolve_homeostasis(HomeostasisRule* homeostasisobject);

		/* Define nested classes to avoid cluttering up the class tree (or file system) */
		class WeightDependentUpdatescalingRule
		{
		private:

			// info for debug output and stuff:
			enum FactoryMethodUsed {Additive,Linear,Guetig2003,Morrison2007,Gilson2011};
			FactoryMethodUsed factory_used;

			WeightDependentUpdatescalingRule() {}	// never let anyone else call the default constructor. Use factory methods instead!
			//WeightDependentUpdatescalingRule(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize);

		protected:
			/** For referring back to the parent connection class. This is set only when adding a WeightDependentUpdateScaling to a STDPwdHomConnection.*/
			STDPwdGrowthConnection* parentConnection;

			/* Constructors: */
			//WeightDependentUpdatescalingRule(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize);
			/** Protected Constructor: Objects of this type should only be produced via the factory methods provided below. */
			WeightDependentUpdatescalingRule(const AurynDouble (WeightDependentUpdatescalingRule::*pFunction_Causal)(const AurynWeight *)const,
										 const AurynDouble (WeightDependentUpdatescalingRule::*pFunction_Anticausal)(const AurynWeight *)const,
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
			const AurynDouble (WeightDependentUpdatescalingRule::*scaleCausal)(const AurynWeight*)const; ///* Function pointer: takes AurynWeight and returns AurynDouble
			const AurynDouble (WeightDependentUpdatescalingRule::*scaleAnticausal)(const AurynWeight*)const; ///* Function pointer: takes AurynWeight and returns AurynDouble

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
			static WeightDependentUpdatescalingRule* makeAdditiveUpdates(     AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize);
			static WeightDependentUpdatescalingRule* makeLinearUpdates(       AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat theAttractorStrengthIndicator, AurynWeight theAttractorLocationIndicator, AurynFloat theMeanSlope=0);
			static WeightDependentUpdatescalingRule* makeGuetig2003Updates(   AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat mu=0.02);
			static WeightDependentUpdatescalingRule* makeMorrison2007Updates( AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat mu_LTP=1, AurynFloat mu_LTD=0.2);
			static WeightDependentUpdatescalingRule* makeGilson2011Updates(   AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat mu_LTP, AurynFloat mu_LTD);

			const AurynDouble dothescale_Causal(const AurynWeight* pWeight)const;
			const AurynDouble dothescale_Anticausal(const AurynWeight* pWeight)const;

			void setMaxWeight(AurynWeight maxweight); ///* To set the max weight before giving this object to the Connection constractor, if needed.
			const string getRuleName()const;

			/** Sets the parent Connection class. This is not done in the constructor to avoid circular dependencies. */
			void setConnection(STDPwdGrowthConnection* pConnection);

		};

		class WeightDependentGrowthRule
		{
		protected:
			STDPwdGrowthConnection* parentConnection; ///* For referring back to the parent connection class. This is set only when adding a GrowthRule to a STDPwdHomConnection for avoiding circular dependencies.

		public:
			/* Class-specific enum definitions: */
			enum TrainednessMethod { Entropy, Kurtosis, SumOfExponentials, SumOfLargeWeights };

			/* Parameters: */
			double strideScale = 0.01;  // default value maxRandomwalkrate=0.01 as in my matlab code (renamed here because the name never fit properly before).
			TrainednessMethod the_trainedness_measure = SumOfLargeWeights;
			bool scaleGrowthStepsByWeight = false;

			AurynFloat updateinterval_trainedness = 1.0f;
			AurynFloat updateinterval_weights = 0.1f;

			/* Internal cache of computed variables (but keep public for debugging for now): */
			//AurynStateVector postsynapticTrainedness;  // object field (no pointer), so that it automatically gets deleted together with this class.
			std::vector<float> postsynapticTrainedness;  // object field (no pointer), so that it automatically gets deleted together with this class.

			/* Constructors & Destructors & Factories: */
			virtual ~WeightDependentGrowthRule() {} ///* empty implementation of virtual destructor so that it can be overridden in subclasses. So no more compiler warnings on this.
			// just the default constructor.
			/** A helper function to easily create growth rules. But you don't need to use it if you don't want to. */
			static WeightDependentGrowthRule* rule_factory(string growth_type, AurynFloat stride, bool scaleByWeight, string trainedness_measure_string);

			/* Implemented (possibly still virtual) functions: */
			/** Sets the parent Connection class. This is not done in the constructor to avoid circular dependencies. */
			void setConnection(STDPwdGrowthConnection* pConnection);

			/** This is the "weight dependent" part of this growth rule. The amount of weight change depends on some function of the set of synaptic weights.
			 * The result of this computation is stored in the field postsynapticTrainedness for later reference. */
			void updateTrainedness();

			/** This actually gets called by other parts of the connection class to grow the weight in each time step. */
			void growWeights();

			/* Abstract virtual functions: */
			virtual AurynFloat get_change() = 0;

			/* Forward declarations; see below */
			class ConstantWDGrowthRule;			// needs to be forward-declared because it inherits from WeightDependentGrowthRule and therefore cannot be defined here yet.
			class GaussRandomWDGrowthRule;		// needs to be forward-declared because it inherits from WeightDependentGrowthRule and therefore cannot be defined here yet.
			class AbsGaussRandomWDGrowthRule;	// needs to be forward-declared because it inherits from WeightDependentGrowthRule and therefore cannot be defined here yet.

		};
	};
	class STDPwdGrowthConnection::WeightDependentGrowthRule::ConstantWDGrowthRule : public STDPwdGrowthConnection::WeightDependentGrowthRule
	{
		virtual AurynFloat get_change();
	};
	class STDPwdGrowthConnection::WeightDependentGrowthRule::GaussRandomWDGrowthRule : public STDPwdGrowthConnection::WeightDependentGrowthRule
	{
		virtual AurynFloat get_change();
	};
	class STDPwdGrowthConnection::WeightDependentGrowthRule::AbsGaussRandomWDGrowthRule : public STDPwdGrowthConnection::WeightDependentGrowthRule
	{
		virtual AurynFloat get_change();
	};

}

#endif /* STDPWDHOMCONNECTION_H_ */
