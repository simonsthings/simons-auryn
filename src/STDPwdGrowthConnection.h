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
		/* Forward declarations of nested classes ; see below */
		class WeightDependentUpdatescalingRule;
		class WeightDependentGrowthRule;

	protected:
		AurynFloat tau_pre;		// use setters to set!
		AurynFloat tau_post;	// use setters to set!
		Trace * tr_pre;
		Trace * tr_post;

		void propagate_forward();
		void propagate_backward();

		AurynWeight get_postspiketrace(NeuronID postspiker);
		AurynWeight get_prespiketrace(NeuronID prespiker);

		//const std::function<AurynDouble(AurynWeight)> &the_weight_dependence;
		WeightDependentUpdatescalingRule* the_weight_dependence = nullptr; // required in constructor.

		//std::vector<WeightDependentGrowthRule*> the_growth_vector; ///< list of homeostasis rules to be applied to the connection in vector's order.
		WeightDependentGrowthRule* the_growth_rule = nullptr;  // init with NULL, set after creation of connection object.

	public:
		bool stdp_active;

		/* Constructors & Destructors & Factories: */
		/** Creates a new connection with a default additive STDP rule. */
		STDPwdGrowthConnection(SpikingGroup *source, NeuronGroup *destination);

		/** Creates a new connection with a given weight dependence rule. */
		STDPwdGrowthConnection(SpikingGroup *source, NeuronGroup *destination, AurynWeight initialweight, AurynWeight maxweight,
									WeightDependentUpdatescalingRule *theWeightDependence, AurynFloat tau_pre=20e-3, AurynFloat tau_post=20e-3,
									AurynFloat sparseness=1.0, TransmitterType transmitter=GLUT);

		/** Deletes the Scaling and Growth rule objects. */
		virtual ~STDPwdGrowthConnection();
		virtual void finalize() override;

		virtual void propagate() override;
		virtual void evolve() override;
		//virtual string get_name() override;

		void setUpdateScalingRule(WeightDependentUpdatescalingRule* pUpdatescalingRule); ///* The update scaling rule can also be changed after object creation, but having it in the constructor ensures that at least one scaling rue is always present.
		void setGrowthRule(WeightDependentGrowthRule* pGrowthRule);		///* If this is NULL at runtime, then just skip growth. No harm done.

		void setTau_pre(AurynFloat tau_pre);
		void setTau_post(AurynFloat tau_post);

	};











	/* Define nested classes to avoid cluttering up the class tree (or file system) */
	class STDPwdGrowthConnection::WeightDependentUpdatescalingRule
	{
	protected:
		/** For referring back to the parent connection class. This is set only when adding a WeightDependentUpdateScaling to a STDPwdHomConnection.*/
		STDPwdGrowthConnection* parentConnection;
		string ruleName; // for quickly returning the rule name for e.g. logging.

		/* fudge (used by every rule): */
		AurynFloat scaleconstant_Causal;
		AurynFloat scaleconstant_Anticausal;

	public:
		/** Sets the parent Connection class. This is not done in the constructor to avoid circular dependencies. */
		void setConnection(STDPwdGrowthConnection* pConnection);
		const virtual string get_name()const;

		/* caller-given parameters (used by every rule): */
		AurynFloat A_plus;			// scale for causal spike order: pre-before-post
		AurynFloat A_minus;			// scale for anticausal spike order: post-before-pre
		AurynFloat update_stepsize; // learning rate

		/* Abstract virtual functions: */
		virtual AurynDouble dothescale_Causal(AurynWeight* pWeight) = 0;		// use more const here?
		virtual AurynDouble dothescale_Anticausal(AurynWeight* pWeight) = 0;

		/* Constructors & Destructors & Factories: */
		virtual ~WeightDependentUpdatescalingRule() {} ///* empty implementation of virtual destructor so that it can be overridden in subclasses. So no more compiler warnings on this.
		WeightDependentUpdatescalingRule(AurynFloat A_plus, AurynFloat A_minus,
										 AurynFloat update_stepsize, string scaling_rule_name);

		/* Forward declarations; see below */
		class AdditiveUpdates;
		class LinearUpdates;
		class Guetig2003Updates;
		class Morrison2007Updates;
	};

	class STDPwdGrowthConnection::WeightDependentUpdatescalingRule::AdditiveUpdates : public STDPwdGrowthConnection::WeightDependentUpdatescalingRule
	{
	public:
		/* rule-specific parameters: */
		/* rule-specific functions: */
		AurynDouble dothescale_Causal(AurynWeight* pWeight) override;
		AurynDouble dothescale_Anticausal(AurynWeight* pWeight) override;
		/** Constructor: */
		AdditiveUpdates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize);
	};
	class STDPwdGrowthConnection::WeightDependentUpdatescalingRule::LinearUpdates : public STDPwdGrowthConnection::WeightDependentUpdatescalingRule
	{
	public:
		/* rule-specific parameters: */
		AurynFloat slope_Causal;
		AurynFloat slope_Anticausal;
		AurynFloat offset_Causal;
		AurynFloat offset_Anticausal;
		/* rule-specific functions: */
		AurynDouble dothescale_Causal(AurynWeight* pWeight) override;
		AurynDouble dothescale_Anticausal(AurynWeight* pWeight) override;
		/** Constructor: */
		LinearUpdates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat theAttractorStrengthIndicator, AurynWeight theAttractorLocationIndicator, AurynFloat theMeanSlope=0);
	};
	class STDPwdGrowthConnection::WeightDependentUpdatescalingRule::Guetig2003Updates : public STDPwdGrowthConnection::WeightDependentUpdatescalingRule
	{
	public:
		/* rule-specific parameters: */
		AurynFloat mu;
		/* rule-specific functions: */
		AurynDouble dothescale_Causal(AurynWeight* pWeight) override;
		AurynDouble dothescale_Anticausal(AurynWeight* pWeight) override;
		/** Constructor: */
		Guetig2003Updates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat mu=0.02);
	};
	class STDPwdGrowthConnection::WeightDependentUpdatescalingRule::Morrison2007Updates : public STDPwdGrowthConnection::WeightDependentUpdatescalingRule
	{
	public:
		/* rule-specific parameters: */
		AurynFloat mu_1;
		AurynFloat mu_2;
		/* rule-specific functions: */
		AurynDouble dothescale_Causal(AurynWeight* pWeight) override;
		AurynDouble dothescale_Anticausal(AurynWeight* pWeight) override;
		/** Constructor: */
		Morrison2007Updates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize, AurynFloat mu_LTP=1, AurynFloat mu_LTD=0.2);
	};












	class STDPwdGrowthConnection::WeightDependentGrowthRule
	{
	protected:
		STDPwdGrowthConnection* parentConnection; ///* For referring back to the parent connection class. This is set only when adding a GrowthRule to a STDPwdHomConnection for avoiding circular dependencies.
		string growthrule_name = "please set me";

	public:
		const string& getGrowthrule_name() const;

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
		/** A helper function to easily create growth rules. But you don't need to use it if you don't want to. Actually, this might need to go away as more options are added... */
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
