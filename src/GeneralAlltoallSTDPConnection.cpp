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

void GeneralAlltoallSTDPConnection::init(AurynFloat eta, AurynFloat tau_pre, AurynFloat tau_post, AurynFloat maxweight, STDPWeightDependence *the_weight_dependence)
{
	if ( dst->get_post_size() == 0 ) return;

	scale_pre = eta; // post-pre
	scale_post = eta; // pre-post

	auryn::logger->parameter("eta",eta);
	auryn::logger->parameter("A",scale_pre);
	auryn::logger->parameter("B",scale_post);

	tr_pre  = src->get_pre_trace(tau_pre);
	tr_post = dst->get_post_trace(tau_post);

	set_min_weight(0.0);
	set_max_weight(maxweight);

	stdp_active = true;
	weight_dependence = new AdditiveWeightDependence(maxweight, 1, 1);
}


void GeneralAlltoallSTDPConnection::finalize() {
	DuplexConnection::finalize();
}

void GeneralAlltoallSTDPConnection::free()
{
}

GeneralAlltoallSTDPConnection::GeneralAlltoallSTDPConnection(SpikingGroup * source, NeuronGroup * destination, TransmitterType transmitter) : DuplexConnection(source, destination, transmitter)
{
}

GeneralAlltoallSTDPConnection::GeneralAlltoallSTDPConnection(SpikingGroup * source, NeuronGroup * destination,
		const char * filename, 
		AurynFloat eta,
		AurynFloat tau_pre,
		AurynFloat tau_post,
		AurynFloat maxweight, 
		TransmitterType transmitter) 
: DuplexConnection(source, 
		destination, 
		filename, 
		transmitter)
{
	init(eta, tau_pre, tau_post, maxweight, nullptr);
}

GeneralAlltoallSTDPConnection::GeneralAlltoallSTDPConnection(SpikingGroup * source, NeuronGroup * destination,
		AurynWeight weight, AurynFloat sparseness, 
		AurynFloat eta, 
		AurynFloat tau_pre,
		AurynFloat tau_post,
		AurynFloat maxweight, 
		TransmitterType transmitter,
		std::string name) 
: DuplexConnection(source, 
		destination, 
		weight, 
		sparseness, 
		transmitter, 
		name)
{
	init(eta, tau_pre, tau_post, maxweight, nullptr);
	if ( name.empty() )
		set_name("GeneralAlltoallSTDPConnection");
}

GeneralAlltoallSTDPConnection::GeneralAlltoallSTDPConnection(SpikingGroup *source, NeuronGroup *destination,
															 AurynWeight initialweight,
															 AurynFloat sparseness,
															 TransmitterType transmitter,
															 string name)
		: DuplexConnection(source, destination, initialweight, sparseness, transmitter, name)
{
	// Todo: put stuff here!
}

GeneralAlltoallSTDPConnection::~GeneralAlltoallSTDPConnection()
{
	if ( dst->get_post_size() > 0 ) 
		free();
	delete weight_dependence;  // or put this into some general auryn cleanup method?
}


AurynWeight GeneralAlltoallSTDPConnection::dw_pre(NeuronID post)
{
	NeuronID translated_spike = dst->global2rank(post); // only to be used for post traces
	AurynDouble dw = scale_pre*tr_post->get(translated_spike);
	return dw;
}

AurynWeight GeneralAlltoallSTDPConnection::dw_post(NeuronID pre)
{
	AurynDouble dw = scale_post*tr_pre->get(pre);
	return dw;
}


void GeneralAlltoallSTDPConnection::propagate_forward()
{
	// loop over all spikes
	for (const NeuronID prespiking : *(src->get_spikes_immediate()) )
	/*
	for (SpikeContainer::const_iterator prespiking = src->get_spikes()->begin() ; // prespiking = pre_spike
			prespiking != src->get_spikes()->end() ; ++prespiking )
	*/
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
			    *weight += dw_pre(*c) * weight_dependence->applyLTDscaling(weight);

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
		/*
		SpikeContainer::const_iterator spikes_end = dst->get_spikes_immediate()->end();
		for (SpikeContainer::const_iterator postspiking = dst->get_spikes_immediate()->begin() ; // postspiking = post_spike
				postspiking != spikes_end ;
				++postspiking )
		*/
		{
			// loop over all presynaptic partners
			for (const NeuronID * c = bkw->get_row_begin(postspiking) ; c != bkw->get_row_end(postspiking) ; ++c ) {

				#ifdef CODE_ACTIVATE_PREFETCHING_INTRINSICS
				// prefetches next memory cells to reduce number of last-level cache misses
				_mm_prefetch((const char *)bkw->get_data(c)+2,  _MM_HINT_NTA);
				#endif

				// computes plasticity update
				AurynWeight * weight = bkw->get_data(c); 
				*weight += dw_post(*c) * weight_dependence->applyLTPscaling(weight);

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

STDPWeightDependence *GeneralAlltoallSTDPConnection::getWeight_dependence() const
{
	return weight_dependence;
}

void GeneralAlltoallSTDPConnection::setWeight_dependence(STDPWeightDependence *weight_dependence)
{
	GeneralAlltoallSTDPConnection::weight_dependence = weight_dependence;
}


