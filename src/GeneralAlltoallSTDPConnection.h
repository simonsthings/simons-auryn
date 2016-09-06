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

#ifndef GENERALIZEDSTDPCONNECTION_H_
#define GENERALIZEDSTDPCONNECTION_H_

#include "auryn/auryn_definitions.h"
#include "auryn/AurynVector.h"
#include "auryn/DuplexConnection.h"
#include "auryn/Trace.h"
#include "auryn/LinearTrace.h"
#include "auryn/SpikeDelay.h"

#include "AdditiveWeightDependence.h"
#include "LinearWeightDependence.h"
#include "Guetig2003WeightDependence.h"
#include "Morrison2007WeightDependence.h"


namespace auryn {


/*! \brief STDP All-to-All Connection with pluggable weight dependence rules.
 *
 * This class implements standard STDP with a double exponential window and optinal
 * offset terms. Window amplitudes, time constants and weight dependence are freely configurable.
 */
class GeneralAlltoallSTDPConnection : public DuplexConnection
{

private:
	void init(AurynFloat eta, AurynFloat tau_pre, AurynFloat tau_post, AurynFloat maxweight, STDPWeightDependence *the_weight_dependence);

protected:

	AurynDouble hom_fudge;

	Trace * tr_pre;
	Trace * tr_post;

	void propagate_forward();
	void propagate_backward();

	AurynWeight dw_pre(NeuronID post);
	AurynWeight dw_post(NeuronID pre);

	STDPWeightDependence* weight_dependence;

public:

	AurynFloat scale_pre; /*!< Amplitude of post-pre part of the STDP window */
	AurynFloat scale_post; /*!< Amplitude of pre-post part of the STDP window */

	bool stdp_active;

	GeneralAlltoallSTDPConnection(SpikingGroup * source, NeuronGroup * destination,
			TransmitterType transmitter=GLUT);

	GeneralAlltoallSTDPConnection(SpikingGroup * source, NeuronGroup * destination,
			const char * filename, 
			AurynFloat eta=1, 
			AurynFloat tau_pre=20e-3,
			AurynFloat tau_post=20e-3,
			AurynFloat maxweight=1. , 
			TransmitterType transmitter=GLUT);

	GeneralAlltoallSTDPConnection(SpikingGroup * source, NeuronGroup * destination,
			AurynWeight weight, AurynFloat sparseness=0.05, 
			AurynFloat eta=1, 
			AurynFloat tau_pre=20e-3,
			AurynFloat tau_post=20e-3,
			AurynFloat maxweight=1. , 
			TransmitterType transmitter=GLUT,
			string name = "GeneralAlltoallSTDPConnection" );

	GeneralAlltoallSTDPConnection(SpikingGroup * source, NeuronGroup * destination,
								  AurynWeight initialweight,
								  AurynFloat sparseness = 1,
								  TransmitterType transmitter=GLUT,
								  string name = "GeneralAlltoallSTDPConnection");

	virtual ~GeneralAlltoallSTDPConnection();
	virtual void finalize();
	void free();

	virtual void propagate();
	virtual void evolve();

	void setWeight_dependence(STDPWeightDependence *weight_dependence);
	STDPWeightDependence *getWeight_dependence() const;
};

}

#endif /* GENERALIZEDSTDPCONNECTION_H_ */
