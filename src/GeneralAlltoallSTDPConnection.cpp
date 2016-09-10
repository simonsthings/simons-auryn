/* 
* Copyright 2014-2016 Friedemann Zenke
*
* This file is part of Auryn, a simulation package for plastic
* spiking neural networks.
* 
* Auryn is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* Auryn is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with Auryn.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "GeneralAlltoallSTDPConnection.h"

using namespace auryn;


GeneralAlltoallSTDPConnection::GeneralAlltoallSTDPConnection(SpikingGroup *source, NeuronGroup *destination)
		: GeneralAlltoallSTDPConnection(source, destination, 0.85, 1.0, new AdditiveWeightDependence())
{
}

GeneralAlltoallSTDPConnection::GeneralAlltoallSTDPConnection(SpikingGroup *source, NeuronGroup *destination, AurynWeight initialweight, AurynWeight maxweight,
															 STDPWeightDependence *theWeightDependence, AurynFloat tau_pre, AurynFloat tau_post, AurynFloat sparseness,
															 TransmitterType transmitter)
		: DuplexConnection(source, destination, initialweight, sparseness, transmitter, "GeneralAlltoallSTDPConnection"), weight_dependence(theWeightDependence)
{
	if ( dst->get_post_size() == 0 ) return;  // not sure how bad this would be if it happens. Taken from class STDPConnection.

	setTau_pre(tau_pre);
	setTau_post(tau_post);

	set_min_weight(0.0);
	set_max_weight(maxweight);

	stdp_active = true;
}

GeneralAlltoallSTDPConnection::~GeneralAlltoallSTDPConnection()
{
	if ( dst->get_post_size() > 0 ) 
		free();
	delete weight_dependence;  // or put this into some general auryn cleanup method?
}

void GeneralAlltoallSTDPConnection::finalize() {
	DuplexConnection::finalize();
}

void GeneralAlltoallSTDPConnection::free()
{
}



AurynWeight GeneralAlltoallSTDPConnection::get_postspikememory(NeuronID postspiker)
{
	NeuronID translated_postspiker = dst->global2rank(postspiker); // only to be used for dst traces
	return tr_post->get(translated_postspiker);
}

AurynWeight GeneralAlltoallSTDPConnection::get_prespikememory(NeuronID prespiker)
{
	return tr_pre->get(prespiker);
}


void GeneralAlltoallSTDPConnection::propagate_forward()
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
			    *weight += get_postspikememory(*c) * weight_dependence->scalePreAfterPost(weight);

				// clips weights
				if ( *weight > get_max_weight() ) *weight = get_max_weight(); 
				else
			    if ( *weight < get_min_weight() ) *weight = get_min_weight();
			}
		}
	}
}

void GeneralAlltoallSTDPConnection::propagate_backward()
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
				*weight += get_prespikememory(*c) * weight_dependence->scalePreBeforePost(weight);

				// clips weights
				if ( *weight > get_max_weight() ) *weight = get_max_weight();
				else
			    if ( *weight < get_min_weight() ) *weight = get_min_weight();
			}
		}
	}
}

void GeneralAlltoallSTDPConnection::propagate()
{
	propagate_forward();
	propagate_backward();
}

void GeneralAlltoallSTDPConnection::evolve()
{
}

void GeneralAlltoallSTDPConnection::setTau_pre(AurynFloat the_tau_pre)
{
	tau_pre = the_tau_pre;
	tr_pre  = src->get_pre_trace(tau_pre);
	auryn::logger->parameter("tau_pre",tau_pre);
}

void GeneralAlltoallSTDPConnection::setTau_post(AurynFloat the_tau_post)
{
	tau_post = the_tau_post;
	tr_post = dst->get_post_trace(tau_post);
	auryn::logger->parameter("tau_post",tau_post);
}

void GeneralAlltoallSTDPConnection::set_max_weight(AurynWeight maximum_weight)
{
	SparseConnection::set_max_weight(maximum_weight);
	weight_dependence->w_max = maximum_weight;
}

string GeneralAlltoallSTDPConnection::get_name()
{
	// Todo: concatenate with the name of the weight dependence:
	return Connection::get_name() + " with " +  weight_dependence->rule_name;
}





