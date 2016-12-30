/*
 * Made by Simon
*/

#include "WDHomeostaticSTDPConnection.h"

using namespace auryn;


WDHomeostaticSTDPConnection::WDHomeostaticSTDPConnection(SpikingGroup *source, NeuronGroup *destination)
		: WDHomeostaticSTDPConnection(source, destination, 0.85, 1.0, WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::makeAdditiveUpdates(0.95,1.0,1.0f/32.0f))
{
}

WDHomeostaticSTDPConnection::WDHomeostaticSTDPConnection(SpikingGroup *source, NeuronGroup *destination, AurynWeight initialweight, AurynWeight maxweight,
														 WeightDependentUpdateScaling *theWeightDependence, AurynFloat tau_pre, AurynFloat tau_post,
														 AurynFloat sparseness, TransmitterType transmitter)
		: DuplexConnection(source, destination, initialweight, sparseness, transmitter, "WDHomeostaticSTDPConnection"), the_weight_dependence(theWeightDependence)
{
	if ( dst->get_post_size() == 0 ) return;  // not sure how bad this would be if it happens. Taken from class STDPConnection.

	setTau_pre(tau_pre);
	setTau_post(tau_post);

	set_min_weight(0.0);
	set_max_weight(maxweight);
	//the_weight_dependence->setMaxWeight(maxweight);

	stdp_active = true;
}

WDHomeostaticSTDPConnection::~WDHomeostaticSTDPConnection()
{
	if ( dst->get_post_size() > 0 ) 
		free();
	delete the_weight_dependence;  // or put this into some general auryn cleanup method?
}

void WDHomeostaticSTDPConnection::finalize() {
	DuplexConnection::finalize();
}

void WDHomeostaticSTDPConnection::free()
{
}




void WDHomeostaticSTDPConnection::setTau_pre(AurynFloat the_tau_pre)
{
	tau_pre = the_tau_pre;
	tr_pre  = src->get_pre_trace(tau_pre);
	auryn::logger->parameter("tau_pre",tau_pre);
}

void WDHomeostaticSTDPConnection::setTau_post(AurynFloat the_tau_post)
{
	tau_post = the_tau_post;
	tr_post = dst->get_post_trace(tau_post);
	auryn::logger->parameter("tau_post",tau_post);
}

/*
void WDHomeostaticSTDPConnection::set_max_weight(AurynWeight maximum_weight)
{
	SparseConnection::set_max_weight(maximum_weight);
	//the_weight_dependence->setMaxWeight(maximum_weight);
}
*/

string WDHomeostaticSTDPConnection::get_name()
{
	// Todo: concatenate with the name of the weight dependence:
	return Connection::get_name() + " with " +  the_weight_dependence->getRuleName() + " weight dependence";
}




AurynWeight WDHomeostaticSTDPConnection::get_postspiketrace(NeuronID postspiker)
{
	NeuronID translated_postspiker = dst->global2rank(postspiker); // only to be used for dst traces
	return tr_post->get(translated_postspiker);
}

AurynWeight WDHomeostaticSTDPConnection::get_prespiketrace(NeuronID prespiker)
{
	return tr_pre->get(prespiker);
}


void WDHomeostaticSTDPConnection::propagate_forward()
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

void WDHomeostaticSTDPConnection::propagate_backward()
{
	if (stdp_active) {
		// loop over all spikes
		for (const NeuronID postspiking : *(dst->get_spikes_immediate()) )
		{
			// loop over all presynaptic partners
			for (const NeuronID * c = bkw->get_row_begin(postspiking) ; c != bkw->get_row_end(postspiking) ; ++c ) {

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

void WDHomeostaticSTDPConnection::propagate()
{
	propagate_forward();
	propagate_backward();
}

void WDHomeostaticSTDPConnection::evolve()
{
	// apply all those non-STDP plasticity things here!
//	for (auto homeostasisobject : homeostasisvector)
	{
		//evolve_homeostasis(homeostasisobject);
	}
}

//void WDHomeostaticSTDPConnection::evolve_homeostasis(HomeostasisRule* homeostasisobject)
//{
//	homeostasisobject->updateWeightsBasedOnWeightDistribution(this);
//}

















const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::scaleCausal_Constant(const AurynWeight *pWeight)const
{
	return scaleconstant_Causal;
}
const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::scaleCausal_Linear(const AurynWeight *pWeight)const
{
	return scaleconstant_Causal * (slope_Causal * (*pWeight) + offset_Causal);
}
const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::scaleCausal_Guetig2003(const AurynWeight *pWeight)const
{
	return scaleconstant_Causal * std::pow(1-*pWeight,mu_1);
}
const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::scaleCausal_Morrison2007(const AurynWeight *pWeight)const
{
	return scaleconstant_Causal * std::pow(*pWeight,mu_1);
}
const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::scaleCausal_Gilson2011(const AurynWeight *pWeight)const
{
	return scaleconstant_Causal * 1;
}


const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::scaleAnticausal_Constant(const AurynWeight *pWeight)const
{
	return scaleconstant_Anticausal;
}
const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::scaleAnticausal_Linear(const AurynWeight *pWeight)const
{
	return scaleconstant_Anticausal * (slope_Anticausal * (*pWeight) + offset_Anticausal);
}
const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::scaleAnticausal_Guetig2003(const AurynWeight *pWeight)const
{
	return scaleconstant_Anticausal * std::pow(*pWeight,mu_2);
}
const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::scaleAnticausal_Morrison2007(const AurynWeight *pWeight)const
{
	return scaleconstant_Anticausal * std::pow(*pWeight,mu_2);
}
const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::scaleAnticausal_Gilson2011(const AurynWeight *pWeight)const
{
	return scaleconstant_Anticausal * 1;
}

/*
WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::WeightDependentUpdateScaling(
		AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize)
		: A_plus(A_plus), A_minus(A_minus),update_stepsize(update_stepsize)
{
	scaleconstant_Causal     = A_plus  * update_stepsize;
	scaleconstant_Anticausal = A_minus * update_stepsize;
}
*/
WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::WeightDependentUpdateScaling(
		const AurynDouble (WeightDependentUpdateScaling::*pFunction_Causal)(const AurynWeight *)const,
		const AurynDouble (WeightDependentUpdateScaling::*pFunction_Anticausal)(const AurynWeight *)const,
		AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize)
		: A_plus(A_plus), A_minus(A_minus),update_stepsize(update_stepsize), scaleCausal(pFunction_Causal), scaleAnticausal(pFunction_Anticausal)
{
	scaleconstant_Causal     = A_plus  * update_stepsize;
	scaleconstant_Anticausal = A_minus * update_stepsize;
}


WDHomeostaticSTDPConnection::WeightDependentUpdateScaling *
WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::makeAdditiveUpdates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize)
{
	// Choose the function pointers so that we don't need to branch during the simulation,
	// and create the weight dependence encapsulation object:
	WeightDependentUpdateScaling* pWDUS = new WeightDependentUpdateScaling(
			&WeightDependentUpdateScaling::scaleCausal_Constant,
			&WeightDependentUpdateScaling::scaleAnticausal_Constant,
			A_plus, A_minus, update_stepsize);

	// assign any special parameters:
	// (nothing to do here)

	// compute any fudge factors:
	// (nothing to do here)

	pWDUS->factory_used = Additive;
	return pWDUS;
}

WDHomeostaticSTDPConnection::WeightDependentUpdateScaling *
WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::makeLinearUpdates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize,
																			 AurynFloat theAttractorStrengthIndicator,
																			 AurynWeight theAttractorLocationIndicator, AurynFloat theMeanSlope)
{
	// Choose the function pointers so that we don't need to branch during the simulation,
	// and create the weight dependence encapsulation object:
	WeightDependentUpdateScaling* pWDUS = new WeightDependentUpdateScaling(
			&WeightDependentUpdateScaling::scaleCausal_Linear,
			&WeightDependentUpdateScaling::scaleAnticausal_Linear,
			A_plus, A_minus, update_stepsize);

	// assign any special parameters:
	//pWDUS->attractorStrengthIndicator = theAttractorStrengthIndicator;
	//pWDUS->attractorLocationIndicator = theAttractorLocationIndicator;

	// compute any fudge factors:
	pWDUS->slope_Causal     = theMeanSlope-theAttractorStrengthIndicator;
	pWDUS->slope_Anticausal = theMeanSlope+theAttractorStrengthIndicator;
	pWDUS->offset_Causal     = 0.5f - pWDUS->slope_Causal     * theAttractorLocationIndicator;
	pWDUS->offset_Anticausal = 0.5f - pWDUS->slope_Anticausal * theAttractorLocationIndicator;

	pWDUS->factory_used = Linear;
	return pWDUS;
}

WDHomeostaticSTDPConnection::WeightDependentUpdateScaling *
WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::makeGuetig2003Updates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize,
																				 AurynFloat mu)
{
	// Choose the function pointers so that we don't need to branch during the simulation,
	// and create the weight dependence encapsulation object:
	WeightDependentUpdateScaling* pWDUS = new WeightDependentUpdateScaling(
			&WeightDependentUpdateScaling::scaleCausal_Guetig2003,
			&WeightDependentUpdateScaling::scaleAnticausal_Guetig2003,
			A_plus, A_minus, update_stepsize);

	// assign any special parameters:
	pWDUS->mu_1 = mu;
	pWDUS->mu_2 = mu;

	// compute any fudge factors:
	// (nothing to do here)

	pWDUS->factory_used = Guetig2003;
	return pWDUS;
}

WDHomeostaticSTDPConnection::WeightDependentUpdateScaling *
WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::makeMorrison2007Updates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize,
																				   AurynFloat mu_LTP, AurynFloat mu_LTD)
{
	// Choose the function pointers so that we don't need to branch during the simulation,
	// and create the weight dependence encapsulation object:
	WeightDependentUpdateScaling* pWDUS = new WeightDependentUpdateScaling(
			&WeightDependentUpdateScaling::scaleCausal_Morrison2007,
			&WeightDependentUpdateScaling::scaleAnticausal_Morrison2007,
			A_plus, A_minus, update_stepsize);

	// assign any special parameters:
	pWDUS->mu_1 = mu_LTP;
	pWDUS->mu_2 = mu_LTD;

	// compute any fudge factors:
	// (nothing to do here)

	pWDUS->factory_used = Morrison2007;
	return pWDUS;
}

WDHomeostaticSTDPConnection::WeightDependentUpdateScaling *
WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::makeGilson2011Updates(AurynFloat A_plus, AurynFloat A_minus, AurynFloat update_stepsize,
																				 AurynFloat mu_LTP, AurynFloat mu_LTD)
{
	WeightDependentUpdateScaling* pWDUS = new WeightDependentUpdateScaling(
			&WeightDependentUpdateScaling::scaleCausal_Gilson2011,
			&WeightDependentUpdateScaling::scaleAnticausal_Gilson2011,
			A_plus, A_minus, update_stepsize);

	// assign any special parameters:
	pWDUS->mu_1 = mu_LTP;
	pWDUS->mu_2 = mu_LTD;
	//TODO: actually implement this according to the paper!

	// compute any fudge factors:
	// (...)

	pWDUS->factory_used = Gilson2011;
	return pWDUS;
}

const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::dothescale_Causal(const AurynWeight* pWeight)const
{
	return (this->*scaleCausal)(pWeight);
}
const AurynDouble WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::dothescale_Anticausal(const AurynWeight* pWeight)const
{
	return (this->*scaleAnticausal)(pWeight);
}

void WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::setMaxWeight(AurynWeight maxweight)
{
	this->compute_fudge(maxweight);
}

const string WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::getRuleName()const
{
	switch (factory_used)
	{
		case Additive:
			return "Additive";
		case Linear:
			return "Linear";
		case Guetig2003:
			return "Guetig2003";
		case Morrison2007:
			return "Morrison2007";
		case Gilson2011:
			return "Gilson2011";
		default:
			throw std::logic_error("This should never happen! Implementation bug?");
	}
}

void WDHomeostaticSTDPConnection::WeightDependentUpdateScaling::compute_fudge(AurynWeight maxweight)
{
	switch (factory_used)
	{
		case Additive:
			// do nothing here.
			break;
		case Linear:
			// TODO: adjust the fudge factors accordingly!
			break;
		case Guetig2003:
			// TODO: Adjust internal max_weight or some fudge.
			break;
		case Morrison2007:
			// TODO: maybe re-scale mu? Or do nothing here?
			break;
		case Gilson2011:
			// TODO implement this.
			break;
		default:
			throw std::logic_error("This should never happen! Implementation bug?");
	}
}




