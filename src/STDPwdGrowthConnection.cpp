/*
 * Made by Simon
*/

#include "STDPwdGrowthConnection.h"

using namespace auryn;


STDPwdGrowthConnection::STDPwdGrowthConnection(SpikingGroup *source, NeuronGroup *destination)
		: STDPwdGrowthConnection(source, destination, 0.85, 1.0, new STDPwdGrowthConnection::WeightDependentUpdatescalingRule::AdditiveUpdates(0.95,1.0,1.0f/32.0f))
{
}

STDPwdGrowthConnection::STDPwdGrowthConnection(SpikingGroup *source, NeuronGroup *destination, AurynWeight initialweight, AurynWeight maxweight,
														 WeightDependentUpdatescalingRule *theWeightDependence, AurynFloat tau_pre, AurynFloat tau_post,
														 AurynFloat sparseness, TransmitterType transmitter)
		: DuplexConnection(source, destination, initialweight, sparseness, transmitter, "STDPwdGrowthConnection")
{
	if ( dst->get_post_size() == 0 ) return;  // not sure how bad this would be if it happens. Taken from class STDPConnection.

	setUpdateScalingRule(theWeightDependence);
	the_weight_dependence->setConnection(this);

	setTau_pre(tau_pre);
	setTau_post(tau_post);

	set_min_weight(0.0);
	set_max_weight(maxweight);

	stdp_active = true;
}

STDPwdGrowthConnection::~STDPwdGrowthConnection()
{
	delete the_weight_dependence;	// or put this into some general auryn cleanup method?
	delete the_growth_rule;			// or put this into some general auryn cleanup method?
}

void STDPwdGrowthConnection::finalize() {
	DuplexConnection::finalize();
}



void STDPwdGrowthConnection::setTau_pre(AurynFloat the_tau_pre)
{
	tau_pre = the_tau_pre;
	tr_pre  = src->get_pre_trace(tau_pre);
	auryn::logger->parameter("tau_pre",tau_pre);
}

void STDPwdGrowthConnection::setTau_post(AurynFloat the_tau_post)
{
	tau_post = the_tau_post;
	tr_post = dst->get_post_trace(tau_post);
	auryn::logger->parameter("tau_post",tau_post);
}


/*
string STDPwdGrowthConnection::get_name()
{
	if (the_growth_rule != NULL)
		return Connection::get_name() + " with " + the_weight_dependence->get_name() + " weight dependence and growth type '"+ the_growth_rule->getGrowthrule_name() +"'.";
	else
		return Connection::get_name() + " with " + the_weight_dependence->get_name() + " weight dependence";
}
*/



AurynWeight STDPwdGrowthConnection::get_postspiketrace(NeuronID postspiker)
{
	NeuronID translated_postspiker = dst->global2rank(postspiker); // only to be used for dst traces
	return tr_post->get(translated_postspiker);
}

AurynWeight STDPwdGrowthConnection::get_prespiketrace(NeuronID prespiker)
{
	return tr_pre->get(prespiker);
}


void STDPwdGrowthConnection::propagate_forward()
{
	// loop over all spikes
	for (const NeuronID prespiking : *(src->get_spikes_immediate()) )
	{
		// loop over all postsynaptic partners
		for (const NeuronID * c = fwd->get_row_begin(prespiking) ; c != fwd->get_row_end(prespiking) ; ++c )
		//for (const NeuronID * c = fwd->get_row_begin(*prespiking) ; c != w->get_row_end(*prespiking) ; ++c )
		{ // c = post index

			// transmit signal to target at postsynaptic neuron
			AurynWeight * weight = fwd->get_data_ptr(c);
			transmit( *c , *weight );

			// handle plasticity
			if ( stdp_active ) {
				// performs weight update
				*weight += get_postspiketrace(*c) * the_weight_dependence->dothescale_Anticausal(weight);

				// clips weights
				if ( *weight > get_max_weight() ) *weight = get_max_weight(); 
				else
			    if ( *weight < get_min_weight() ) *weight = get_min_weight();
			}
		}
	}
}

void STDPwdGrowthConnection::propagate_backward()
{
	if (stdp_active) {
		// loop over all spikes
		for (const NeuronID postspiking : *(dst->get_spikes_immediate()) )
		{
			// loop over all presynaptic partners
			for (const NeuronID * c = bkw->get_row_begin(postspiking) ; c != bkw->get_row_end(postspiking) ; ++c )
			{

				#ifdef CODE_ACTIVATE_PREFETCHING_INTRINSICS
				// prefetches next memory cells to reduce number of last-level cache misses
				_mm_prefetch((const char *)bkw->get_data(c)+2,  _MM_HINT_NTA);
				#endif

				// computes plasticity update
				AurynWeight * weight = bkw->get_data(c); 
				*weight += get_prespiketrace(*c) * the_weight_dependence->dothescale_Causal(weight);

				// clips weights
				if ( *weight > get_max_weight() ) *weight = get_max_weight();
				else
			    if ( *weight < get_min_weight() ) *weight = get_min_weight();
			}
		}
	}
}

void STDPwdGrowthConnection::propagate()
{
	propagate_forward();
	propagate_backward();
}

void STDPwdGrowthConnection::evolve()
{
	// apply all those non-STDP plasticity things here!
	WeightDependentGrowthRule* a_homeostasis_rule = the_growth_rule;
	//for (auto a_homeostasis_rule : homeostasisvector)
	if (the_growth_rule != nullptr)
	{
		//evolve_homeostasis(a_homeostasis_rule);
		a_homeostasis_rule->growWeights();
	}
}

void STDPwdGrowthConnection::setUpdateScalingRule(WeightDependentUpdatescalingRule* pUpdatescalingRule)
{
	if (pUpdatescalingRule != NULL)
	{
		// delete any old rule object:
		if (the_weight_dependence != NULL)
			delete the_weight_dependence;

		the_weight_dependence = pUpdatescalingRule;
		the_weight_dependence->setConnection(this);

		string namestring = "STDPwdGrowthConnection with " + the_weight_dependence->get_name() + " weight dependence.";
		if (the_growth_rule != NULL)
			namestring = "STDPwdGrowthConnection with " + the_weight_dependence->get_name() + " weight dependence and growth type '"+ the_growth_rule->getGrowthrule_name() +"'.";
		set_name(namestring);
	}
	else
		throw std::invalid_argument( "The update scaling rule must never be NULL. Please change your implementation." );
}

void STDPwdGrowthConnection::setGrowthRule(WeightDependentGrowthRule* pGrowthRule)
{
	// delete any old rule object:
	if (the_growth_rule != NULL)
		delete the_growth_rule;

	the_growth_rule = pGrowthRule;
	string namestring = "STDPwdGrowthConnection with " + the_weight_dependence->get_name() + " weight dependence.";
	if (the_growth_rule != NULL)
	{
		the_growth_rule->setConnection(this);
		namestring = "STDPwdGrowthConnection with " + the_weight_dependence->get_name() + " weight dependence and growth type '"+ the_growth_rule->getGrowthrule_name() +"'.";
	}
	set_name(namestring);

}












STDPwdGrowthConnection::WeightDependentUpdatescalingRule::WeightDependentUpdatescalingRule(AurynFloat A_plus, AurynFloat A_minus,
																						   AurynFloat update_stepsize, string scaling_rule_name)
		: ruleName(scaling_rule_name),A_plus(A_plus), A_minus(A_minus),update_stepsize(update_stepsize)
{
	scaleconstant_Causal     = A_plus  * update_stepsize;
	scaleconstant_Anticausal = A_minus * update_stepsize;
}

const string STDPwdGrowthConnection::WeightDependentUpdatescalingRule::get_name()const
{
	return ruleName;
}

void STDPwdGrowthConnection::WeightDependentUpdatescalingRule::setConnection(STDPwdGrowthConnection* pConnection)
{
	parentConnection = pConnection;
}

AurynDouble STDPwdGrowthConnection::WeightDependentUpdatescalingRule::AdditiveUpdates::dothescale_Causal(AurynWeight* pWeight)
{
	return scaleconstant_Causal;
}
AurynDouble STDPwdGrowthConnection::WeightDependentUpdatescalingRule::LinearUpdates::dothescale_Causal(AurynWeight* pWeight)
{
	return scaleconstant_Causal * (slope_Causal * (*pWeight) + offset_Causal);
}
AurynDouble STDPwdGrowthConnection::WeightDependentUpdatescalingRule::Guetig2003Updates::dothescale_Causal(AurynWeight* pWeight)
{
	return scaleconstant_Causal * std::pow(1-*pWeight,mu);
}
AurynDouble STDPwdGrowthConnection::WeightDependentUpdatescalingRule::Morrison2007Updates::dothescale_Causal(AurynWeight* pWeight)
{
	return scaleconstant_Causal * std::pow(*pWeight,mu_1);
}

AurynDouble STDPwdGrowthConnection::WeightDependentUpdatescalingRule::AdditiveUpdates::dothescale_Anticausal(AurynWeight* pWeight)
{
	return scaleconstant_Anticausal;
}
AurynDouble STDPwdGrowthConnection::WeightDependentUpdatescalingRule::LinearUpdates::dothescale_Anticausal(AurynWeight* pWeight)
{
	return scaleconstant_Anticausal * (slope_Anticausal * (*pWeight) + offset_Anticausal);
}
AurynDouble STDPwdGrowthConnection::WeightDependentUpdatescalingRule::Guetig2003Updates::dothescale_Anticausal(AurynWeight* pWeight)
{
	return scaleconstant_Anticausal * std::pow(*pWeight,mu);
}
AurynDouble STDPwdGrowthConnection::WeightDependentUpdatescalingRule::Morrison2007Updates::dothescale_Anticausal(AurynWeight* pWeight)
{
	return scaleconstant_Anticausal * std::pow(*pWeight,mu_2);
}


STDPwdGrowthConnection::WeightDependentUpdatescalingRule::AdditiveUpdates::AdditiveUpdates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize)
		: WeightDependentUpdatescalingRule(A_plus,A_minus,update_stepsize,"Additive")
{
	// no fudge to compute here for additive STDP rules.
}
STDPwdGrowthConnection::WeightDependentUpdatescalingRule::LinearUpdates::LinearUpdates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize,
		AurynFloat theAttractorStrengthIndicator, AurynWeight theAttractorLocationIndicator, AurynFloat theMeanSlope)
		: WeightDependentUpdatescalingRule(A_plus,A_minus,update_stepsize,"Linear")
{
	slope_Causal     = theMeanSlope-theAttractorStrengthIndicator;
	slope_Anticausal = theMeanSlope+theAttractorStrengthIndicator;
	offset_Causal     = 0.5f - slope_Causal     * theAttractorLocationIndicator;
	offset_Anticausal = 0.5f - slope_Anticausal * theAttractorLocationIndicator;

}
STDPwdGrowthConnection::WeightDependentUpdatescalingRule::Guetig2003Updates::Guetig2003Updates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize,
	   AurynFloat mu)
		: WeightDependentUpdatescalingRule(A_plus,A_minus,update_stepsize,"Guetig2003"),
		  mu(mu)
{
	// no fudge to compute here for Guetig2003 STDP rules.
}
STDPwdGrowthConnection::WeightDependentUpdatescalingRule::Morrison2007Updates::Morrison2007Updates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize,
	   AurynFloat mu_LTP, AurynFloat mu_LTD)
		: WeightDependentUpdatescalingRule(A_plus,A_minus,update_stepsize,"Morrison2007"),
		  mu_1(mu_LTP),mu_2(mu_LTD)
{
	// no fudge to compute here for Morrison2007 STDP rules.
}


STDPwdGrowthConnection::WeightDependentGrowthRule*
STDPwdGrowthConnection::WeightDependentGrowthRule::rule_factory(string growth_type, AurynFloat stride, bool scaleByWeight, string trainedness_measure_string)
{
	// options: None, ConstantGrowth, RandomGrowth, RandomJitter, RandomShrinkage, ConstantShrinkage
	if (growth_type == "None")
		return nullptr;
	else
	{
		STDPwdGrowthConnection::WeightDependentGrowthRule* the_growth_rule;
		if (growth_type == "ConstantGrowth")
			the_growth_rule = new STDPwdGrowthConnection::WeightDependentGrowthRule::ConstantWDGrowthRule();
		else if (growth_type == "RandomGrowth")
			the_growth_rule = new STDPwdGrowthConnection::WeightDependentGrowthRule::AbsGaussRandomWDGrowthRule();
		else if (growth_type == "RandomJitter")
			the_growth_rule = new STDPwdGrowthConnection::WeightDependentGrowthRule::GaussRandomWDGrowthRule();
		else if (growth_type == "RandomShrinkage")
			throw std::invalid_argument("not implemented");
		else if (growth_type == "ConstantShrinkage")
			throw std::invalid_argument("not implemented");
		else
			throw std::invalid_argument( "Unknown type of growth rule requested. Typo?" );

		the_growth_rule->growthrule_name = growth_type;
		the_growth_rule->strideScale = stride;
		the_growth_rule->scaleGrowthStepsByWeight = scaleByWeight;

		// Define trainedness measure from the given string:
		// options: Entropy, Kurtosis, SumOfExponentials, SumOfLargeWeights
		if (trainedness_measure_string == "Entropy")
			the_growth_rule->the_trainedness_measure = TrainednessMethod::Entropy;
		else if (trainedness_measure_string == "Kurtosis")
			the_growth_rule->the_trainedness_measure = TrainednessMethod::Kurtosis;
		else if (trainedness_measure_string == "SumOfExponentials")
			the_growth_rule->the_trainedness_measure = TrainednessMethod::SumOfExponentials;
		else if (trainedness_measure_string == "SumOfLargeWeights")
			the_growth_rule->the_trainedness_measure = TrainednessMethod::SumOfLargeWeights;
		else
			throw std::invalid_argument( "Unknown type of growth rule trainedness measure requested! Typo?)" );

		return the_growth_rule;
	}
}


void STDPwdGrowthConnection::WeightDependentGrowthRule::updateTrainedness()
{
	NeuronID N_pre  = parentConnection->src->get_size();
	NeuronID N_post = parentConnection->dst->get_size();

	/* Switch over the different possible Traindness Measures as an enum for now. Perhaps get rid of the If statements if this really gives large performence gains?*/
	switch (the_trainedness_measure)
	{
		//case Entropy:
		//case Kurtosis:
		//case SumOfExponentials:
		case SumOfLargeWeights:
		{
			float strongWeightsSumGoal = 30.0f;
			double strongVal = 0.75; // * maxWeight
			// compute trainedness for each postsynaptic partner (assuming dense connection matrix):
			for (NeuronID postID = 0 ; postID < N_post ; postID++ )
			{
				unsigned int strongSum = 0;  // counter for the number of strong connections for this postsynaptic neuron.
				for (NeuronID preID = 0 ; preID < N_pre ; preID++ )
				{
					const AurynWeight& weight = parentConnection->w->get(preID, postID);
					if (weight >= strongVal) strongSum++;
				}
				float proportionOfStrongConnections = strongSum / strongWeightsSumGoal;
				if (proportionOfStrongConnections > 1)
					//this->postsynapticTrainedness.set(postID,1);
					this->postsynapticTrainedness[postID] = 1;
				else
					//this->postsynapticTrainedness.set(postID,proportionOfStrongConnections);
					this->postsynapticTrainedness[postID] = proportionOfStrongConnections;
			}
		}
			break;
		default:
			throw std::logic_error("This should never happen! Implementation missing or bug?");
	}
}

void STDPwdGrowthConnection::WeightDependentGrowthRule::setConnection(STDPwdGrowthConnection* pConnection)
{
	parentConnection = pConnection;

	// ensure that the vector has the correct size (this will likely be constant-runtime in all but the very first call, as no more resizing will be necessary)
	NeuronID N_post = parentConnection->dst->get_size();
	postsynapticTrainedness.resize(N_post);
}


void STDPwdGrowthConnection::WeightDependentGrowthRule::growWeights()
{
	// only update trainedness if it is time to do so:
	if ( auryn::sys->get_clock() % AurynTime(updateinterval_trainedness/auryn_timestep ) == 0 ) // potential for optimisation by pre-casting and pre-computing, but would require own constructor.
	{
		/* This is the "weight dependent" part of this growth rule. Results stored in field postsynapticTrainedness for reference below. */
		updateTrainedness();
	}

	// only update the weights if it is time to do so:
	if ( auryn::sys->get_clock() % AurynTime(updateinterval_weights/auryn_timestep ) == 0 ) // potential for optimisation by pre-casting and pre-computing, but would require own constructor.
	{
		//NeuronID N_pre = parentConnection->src->get_size();
		NeuronID N_post = parentConnection->dst->get_size();

		// vector operation in matlab: may need to be put into a loop here? Or apply the new vectorized functions of Auryn?
		for (NeuronID postID = 0 ; postID < N_post ; postID++ )
		{
			//double distribution_dependent_weight_change = strideScale * (1 - postsynapticTrainedness.get(postID) );
			double distribution_dependent_weight_change = strideScale * (1 - postsynapticTrainedness[postID] );

			// loop over all presynaptic partners
			for (const NeuronID * c = parentConnection->bkw->get_row_begin(postID) ; c != parentConnection->bkw->get_row_end(postID) ; ++c )
			{
				#ifdef CODE_ACTIVATE_PREFETCHING_INTRINSICS
				// prefetches next memory cells to reduce number of last-level cache misses
				_mm_prefetch((const char *)parentConnection->bkw->get_data(c)+2,  _MM_HINT_NTA);
				#endif

				// computes weight update
				AurynWeight* pWeight = parentConnection->bkw->get_data(c);
				if (scaleGrowthStepsByWeight)
					// just scale by the mean of the STDP weight dependence for this weight. (We could probably also distinguish between weight increases and decreases, but... who cares right now? And would that even make sense?)
					*pWeight += distribution_dependent_weight_change * get_change() * (parentConnection->the_weight_dependence->dothescale_Causal(pWeight) + parentConnection->the_weight_dependence->dothescale_Anticausal(pWeight)) / 2.0f;
				else
					*pWeight += distribution_dependent_weight_change * get_change();


				// clips weights
				if      ( *pWeight > parentConnection->get_max_weight() ) *pWeight = parentConnection->get_max_weight();
				else if ( *pWeight < parentConnection->get_min_weight() ) *pWeight = parentConnection->get_min_weight();
			}
		}
	}
}

const string& STDPwdGrowthConnection::WeightDependentGrowthRule::getGrowthrule_name() const
{
	return growthrule_name;
}


AurynFloat STDPwdGrowthConnection::WeightDependentGrowthRule::ConstantWDGrowthRule::get_change()
{
	return 1;
}

AurynFloat STDPwdGrowthConnection::WeightDependentGrowthRule::GaussRandomWDGrowthRule::get_change()
{
//	return randn(1);
	throw std::logic_error("This has not been implemented yet.");
	/* Matlab code:	  randomwalkeffect = this.maxRandomwalkrate .* bsxfun(@times, (1-this.postsynapticTrainedness), randn(size(this.weights)) );	*/
}
AurynFloat STDPwdGrowthConnection::WeightDependentGrowthRule::AbsGaussRandomWDGrowthRule::get_change()
{
//	return abs(randn(1));
	throw std::logic_error("This has not been implemented yet.");
}

