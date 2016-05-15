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

#ifndef IZHI2003GROUP_H_
#define IZHI2003GROUP_H_

#include "auryn_definitions.h"
#include "NeuronGroup.h"
#include "System.h"

enum Izhikevich_DefaultParametersets {default_2003, fitted_2003, fig2004_A, fig2004_B, fig2004_C, fig2004_D, fig2004_E, fig2004_F, fig2004_G
	, fig2004_H, fig2004_I, fig2004_J, fig2004_K, fig2004_L, fig2004_M, fig2004_N, fig2004_O, fig2004_P, fig2004_Q, fig2004_R, fig2004_S, fig2004_T
	, book2007_RS};

/*! \brief Implements the Izhikevich neuron model as published in (2003,2004,2007).
 *
 * This model hopefully resolves some inconsistencies between the different published versions.
 */
class Izhikevich2003Group : public NeuronGroup
{
private:
	auryn_vector_float * t_leak;
	auryn_vector_float * t_exc;
	auryn_vector_float * t_inh;
    AurynFloat * handle_mem;
    AurynFloat * handle_u;
    AurynFloat * handle_u_temp;
    AurynFloat * handle_e_input;
	AurynFloat scale_mem, scale_u;

	//AurynFloat scale_ampa,scale_gaba, scale_thr;
	//AurynFloat tau_ampa,tau_gaba,tau_nmda;
	//AurynFloat A_ampa,A_nmda;
	//AurynFloat e_rest,e_rev,thr_rest,tau_mem,tau_thr,dthr;

	// SI unit converters:
	AurynFloat mV = 1e-3;  // needed to generalise Izhikevich model to SI units and other timestep durations.
	AurynFloat ms = 1e-3;  // needed to generalise Izhikevich model to SI units and other timestep durations.
	AurynFloat pF = 1e-12; // needed to generalise Izhikevich model to SI units and other timestep durations.

	// izhikevich parameters:
	AurynDouble mem_coeffs[3]; // used to pre-compute some stuff and make the 2003 and 2007 versions compatible.
	AurynDouble u_coeffs[2];       // needed to make the 2003, 2004 and 2007 versions of Izhikevich code compatible.

	AurynFloat izhi_a,izhi_b,izhi_c,izhi_d;
	AurynFloat V_reset, V_init, U_init;
	AurynFloat V_peak; // spike top cutoff at +30mV. This is not a LIF spiek threshold, as the izhikevich neuron produces its own spike upstrokes.
	AurynFloat V_min; // used to stabilize izhikevich neurons! Otherwise, strong input pulses might push it into extreme firing!

	AurynFloat C; // capacitance C as in Izhikevich book p. 273 (Chapter 8.1.4)
	AurynFloat k; // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
	AurynFloat v_rest; // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
	AurynFloat v_thres; // soft threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)

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
	void calculate_Izhikevich_coefficients(AurynDouble k, AurynDouble v_rest, AurynDouble v_thres, AurynDouble b);
	void integrate_membrane();
	void integrate_membrane_debug();
	void integrate_membrane_debug_two();
	void check_peaks();

public:
	bool use2003integration = false; ///< a switch to exactly reproduce the behaviour of the 2003 matlab code. That code was numerically not quite correct, but was used a lot due to easy sharing. Not yet implemented.

	/*! Default constructor.
	 *
	 * @param size the size of the group.  @param load a load specifier that
	 * can in principle help to assign the right amount of neurons to a single
	 * node. However this feature is experimental.  
	 * @param total the total number of neurons that are going to be assigned
	 * in this simulation. This is also part of the exerimental load balancing
	 * mechanism.
	 */
	Izhikevich2003Group( NeuronID size, AurynFloat load = 1.0, NeuronID total = 0 );
	virtual ~Izhikevich2003Group();
	void clear();

	/**
	 * For reproducing the overview image of Izhikevich (2004) and Figure 8.8 (p.277) of the book, we provide the default example settings
	 * for this neuron model.
	 */
	void select_example_parameter_set(Izhikevich_DefaultParametersets paramsetchoice);

	/*! Internally used evolve function. Called by System. */
	virtual void evolve();

	AurynState get_t_exc(NeuronID i);
	AurynState get_t_inh(NeuronID i);
	AurynState get_v_temp(NeuronID i);
	AurynState get_u_temp(NeuronID i);
	AurynState get_u(NeuronID i);
	AurynState get_tempMemState(int i);


};

#endif /*IZHI2003GROUP_H_*/

