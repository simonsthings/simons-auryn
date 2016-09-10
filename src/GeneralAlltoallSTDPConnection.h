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

	AurynWeight get_postspikememory(NeuronID postspiker);
	AurynWeight get_prespikememory(NeuronID prespiker);

	STDPWeightDependence* weight_dependence;

public:
	bool stdp_active;

	/** Creates a new connection with a default additive STDP rule. */
	GeneralAlltoallSTDPConnection(SpikingGroup *source, NeuronGroup *destination);

	/** Creates a new connection with a given weight dependence rule. */
	GeneralAlltoallSTDPConnection(SpikingGroup *source, NeuronGroup *destination, AurynWeight initialweight, AurynWeight maxweight,
								  STDPWeightDependence *theWeightDependence, AurynFloat tau_pre=20e-3, AurynFloat tau_post=20e-3, AurynFloat sparseness=1.0,
								  TransmitterType transmitter=GLUT);

	virtual ~GeneralAlltoallSTDPConnection();
	virtual void finalize() override;
	void free();

	virtual void propagate() override;
	virtual void evolve() override;

	virtual string get_name() override;

	virtual void set_max_weight(AurynWeight maximum_weight) override;

	void setTau_pre(AurynFloat tau_pre);
	void setTau_post(AurynFloat tau_post);

};

}

#endif /* GENERALIZEDSTDPCONNECTION_H_ */
