/* 
* Copyright 2014-2015 Friedemann Zenke, Izhikevich Model implemented by Simon
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
*
* If you are using Auryn or parts of it for your work please cite:
* Zenke, F. and Gerstner, W., 2014. Limits to high-speed simulations 
* of spiking neural networks using general-purpose computers. 
* Front Neuroinform 8, 76. doi: 10.3389/fninf.2014.00076
*/

#ifndef IZHIGROUP_H_
#define IZHIGROUP_H_

#include "auryn_definitions.h"
#include "NeuronGroup.h"
#include "System.h"


/*! \brief Implements the standard integrate and file model used in Auryn.
 *
 * This is the implementation of a standard leaky integrate and fire model
 * with relative refractoriness and conductance based synapses. The model
 * has two time constants for excitatory conductances (AMPA,NMDA). AMPA 
 * conductances are modelled as exponential decays. NMDA is implemented to 
 * have a double exponental timecourse, by low-pass filtering over the AMPA
 * equation. The amplitude between the individual contributions can be 
 * ajusted via set_ampa_nmda_ratio. The voltage dependence of NMDA is 
 * ignored in this model.
 */
class IzhikevichGroup : public NeuronGroup
{
private:
	auryn_vector_float * t_leak;
	auryn_vector_float * t_exc;
	auryn_vector_float * t_inh;
	AurynFloat scale_ampa,scale_gaba, scale_thr;
	AurynFloat tau_ampa,tau_gaba,tau_nmda;
	AurynFloat A_ampa,A_nmda;
	AurynFloat e_rest,e_rev,thr_rest,tau_mem,tau_thr,dthr;

	// izhikevich parameters:
	AurynFloat mV,ms; // needed to generalise Izhikevich model to SI units and other timestep durations.
	AurynFloat a,b,c,d;
	AurynFloat V_peak; // spike top cutoff at +30mV. This is not a LIF spiek threshold, as the izhikevich neuron produces its own spike upstrokes.
	AurynFloat V_min; // used to stabilize izhikevich neurons! Otherwise, strong input pulses might push it into extreme firing!

	AurynFloat C; // capacitance C as in Izhikevich book p. 273 (Chapter 8.1.4)
	AurynFloat k; // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
	AurynFloat v_rest; // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
	AurynFloat v_thres; // soft threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)

	bool doProperChannelStuff;
	bool consistent_integration;
	bool use_recovery; // one-dimensional model as mentioned in Izhikevich (2004) "Choosing models" paper.
	auryn_vector_float * v_temp;
	auryn_vector_float * v_temp2; // using second temp var for computations for now. May be optimised in the future.
	auryn_vector_float * u_temp;
	auryn_vector_float * u;

	auryn_vector_float * inputCurrents;
	auryn_vector_float * backgroundCurrents;

	// debugging only:
	float tempMemStates[20];

	void init();
	void free();
	void calculate_scale_constants();
	void integrate_linear_nmda_synapses();
	void integrate_membrane();
	void integrate_membrane_debug();
	void check_peaks();
	//void integrate_membrane_IFNeuron();
public:
	AurynFloat e_reset;
	/*! Default constructor.
	 *
	 * @param size the size of the group.  @param load a load specifier that
	 * can in principle help to assign the right amount of neurons to a single
	 * node. However this feature is experimental.  
	 * @param total the total number of neurons that are going to be assigned
	 * in this simulation. This is also part of the exerimental load balancing
	 * mechanism.
	 */
	IzhikevichGroup( NeuronID size, AurynFloat load = 1.0, NeuronID total = 0 );
	virtual ~IzhikevichGroup();
	/*! Sets the membrane time constant */
	void set_tau_mem(AurynFloat taum);
	/*! Returns the membrane time constant */
	AurynFloat get_tau_mem();

	/*! Sets the exponential decay time constant of the AMPA conductance (default=5ms). */
	void set_tau_ampa(AurynFloat tau);

	/*! Returns the exponential decay time constant of the AMPA conductance. */
	AurynFloat get_tau_ampa();

	/*! Sets the exponential decay time constant of the GABA conductance (default=10ms). */
	void set_tau_gaba(AurynFloat tau);

	/*! Returns the exponential decay time constant of the GABA conductance. */
	AurynFloat get_tau_gaba();

	/*! Sets the exponential decay time constant of the NMDA conductance (default=100ms).
	 * The rise is governed by tau_ampa if tau_ampa << tau_nmda. */
	void set_tau_nmda(AurynFloat tau);

	/*! Returns the exponential decay time constant of the NMDA conductance.
	 * The rise is governed by tau_ampa if tau_ampa << tau_nmda. */
	AurynFloat get_tau_nmda();
	/*! Set ratio between ampa/nmda contribution to excitatory conductance. */
	void set_ampa_nmda_ratio(AurynFloat ratio);
	void clear();
	/*! Internally used evolve function. Called by System. */
	virtual void evolve();


	AurynState get_t_exc(NeuronID i);
	AurynState get_t_inh(NeuronID i);
	AurynState get_v_temp(NeuronID i);
	AurynState get_u_temp(NeuronID i);
	AurynState get_u(NeuronID i);
	AurynState get_tempMemState(int i);


};

#endif /*IZHIGROUP_H_*/

