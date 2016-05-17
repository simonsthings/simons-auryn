/* 
* Copyright 2014-2015 Friedemann Zenke
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

#include <Izhikevich2003Group.h>


Izhikevich2003Group::Izhikevich2003Group( NeuronID size, AurynFloat load, NeuronID total ) : NeuronGroup(size,load,total)
{
	sys->register_spiking_group(this);
	if ( evolve_locally() ) init();
}

void Izhikevich2003Group::calculate_scale_constants()
{
}

void Izhikevich2003Group::init()
{
	// BEGIN old stuff
	calculate_scale_constants();
	t_leak = get_state_vector("t_leak");
	t_exc =  get_state_vector("t_exc");
	t_inh = get_state_vector("t_inh");
	// END old stuff
	
	// new stuff!

	//v = get_state_vector("v");  // using mem instead!
	u = get_state_vector("u");
	v_temp = get_state_vector("v_temp");
	v_temp2 = get_state_vector("v_temp2");
	u_temp = get_state_vector("u_temp");
	inputCurrents =  get_state_vector("inputCurrents");
	backgroundCurrents = get_state_vector("backgroundCurrents");
	consistent_integration = true;
	use_recovery = false;

	select_example_parameter_set(Izhikevich_DefaultParametersets::default_2003);
	cout << "The Izhikevich coefficients are now: " << endl;
	cout << " a2 = " << mem_coeffs[2];
	cout << ", a1 = " << mem_coeffs[1];
	cout << ", a0 = " << mem_coeffs[0] << endl;
	cout << " U:     a1 = " << u_coeffs[1];
	cout << ", a0 = " << u_coeffs[0] << endl;

	/** Begin set up Handles **/
    handle_mem = auryn_vector_float_ptr ( mem , 0 );
    handle_u = auryn_vector_float_ptr ( u , 0 );
    handle_u_temp = auryn_vector_float_ptr ( u_temp , 0 );
    handle_e_input = auryn_vector_float_ptr ( t_exc , 0 );
    /** End set up Handles **/

    // scale constants to keep time and other units SI-conform:
	scale_mem = dt/ms /(C/pF) * mV;
	scale_u = dt/ms *mV;

	clear();
}
void Izhikevich2003Group::clear()
{
	clear_spikes();

	//auryn_vector_float_set (mem, i, e_rest);
	auryn_vector_float_set_all(mem, V_init);
	auryn_vector_float_set_all (thr, 0.);
	auryn_vector_float_set_all (g_ampa, 0.);
	auryn_vector_float_set_all (g_gaba, 0.);
	auryn_vector_float_set_all (g_nmda, 0.);

	auryn_vector_float_set_all (u, U_init);
	auryn_vector_float_set_all (inputCurrents, 0.);
	auryn_vector_float_set_all (backgroundCurrents, 0.);
}

void Izhikevich2003Group::select_example_parameter_set(Izhikevich_DefaultParametersets paramsetchoice)
{
	//AurynDouble k,v_rest,v_thres;
	V_min = -150; // -150 Volts. This is just a safeguard. TODO: Output warning if this is reached!

	switch (paramsetchoice)
	{
	case Izhikevich_DefaultParametersets::default_2003:
		izhi_a = 0.02;
		izhi_b = 0.2;
		izhi_c = -65; // mV, but don't use SI for now.
		izhi_d = 6;
		// To allow (0.04v^2 + 5v + 140) as in Izhikevich (2003) paper, use these params:
		C = 1;
		k = 0.04;                   // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_thres = -42.3444;//*mV;  // soft threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_rest  = -82.6556;//*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)

		V_peak = 30*mV;
		V_reset = izhi_c*mV;
		V_init = -70*mV;
		U_init = izhi_b*V_init/mV;
		mem_coeffs[0] = 140;
		mem_coeffs[1] = 5 /mV;
		mem_coeffs[2] = 0.04 /mV/mV;
		u_coeffs[0] = 0;
		u_coeffs[1] = izhi_b /mV;
		break;
	case Izhikevich_DefaultParametersets::fitted_2003:
		izhi_a = 0.02;
		izhi_b = 0.2;
		izhi_c = -65;
		izhi_d = 6;
		// To allow (0.04v^2 + 5v + 140) as in Izhikevich (2003) paper, use these params:
		C = 1;
		k = 0.04;                   // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_thres = -42.3444;//*mV;  // soft threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_rest  = -82.6556;//*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)

		V_peak = 30*mV;
		V_reset = izhi_c*mV;
		V_init = -70*mV;
		U_init = izhi_b*V_init/mV;
		calculate_Izhikevich_coefficients(k,v_rest,v_thres,izhi_b);
		u_coeffs[0] = 0;
		break;
	case Izhikevich_DefaultParametersets::fig2004_A:
		// To allow RS neuron of p. 274 of Izhikevich book, use these params:
		izhi_a = 0.03;
		izhi_b = -2;
		izhi_c = -50;
		izhi_d = 100;
		C = 100;  // pF
		k = 0.7;                // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_rest  = -60;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_thres = -40;  // instantaneous threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		// I_syn = 70; // pA
		V_peak = 35*mV;
		V_reset = izhi_c*mV;
		V_init = v_rest*mV;
		U_init = 0;
		calculate_Izhikevich_coefficients(k,v_rest,v_thres,izhi_b);
		break;
	case Izhikevich_DefaultParametersets::book2007_RS:
		// To allow RS neuron of p. 274 of Izhikevich book, use these params:
		izhi_a = 0.03;
		izhi_b = -2;
		izhi_c = -50;
		izhi_d = 100;
		C = 100;  // pF
		k = 0.7;                // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_rest  = -60;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_thres = -40;  // instantaneous threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		// I_syn = 70; // pA
		V_peak = 35*mV;
		V_reset = izhi_c*mV;
		V_init = v_rest*mV;
		U_init = 0;
		calculate_Izhikevich_coefficients(k,v_rest,v_thres,izhi_b);
		break;
	default:
		throw std::logic_error("Unhandled switch case!");
	}

}
void Izhikevich2003Group::calculate_Izhikevich_coefficients(AurynDouble k, AurynDouble v_rest, AurynDouble v_thres, AurynDouble b)
{
	mem_coeffs[0] = k * v_rest * v_thres;
	mem_coeffs[1] = k * -(v_rest + v_thres)  /mV;
	mem_coeffs[2] = k  /mV/mV;

	u_coeffs[0] = b*-v_rest;
	u_coeffs[1] = b   /mV;
}
void Izhikevich2003Group::free()
{
}
Izhikevich2003Group::~Izhikevich2003Group()
{
	if ( evolve_locally() ) free();
}
void Izhikevich2003Group::integrate_membrane_debug_two()
{
	AurynFloat* handle_inputcurrent = auryn_vector_float_ptr ( t_exc , 0 );

	if (use_recovery)
	{
		// will refer to this below via handle_u_temp:
		auryn_vector_float_copy(u,u_temp);
	}

    // TODO we should vectorize this code and use some fast SSE
    // library such as http://gruntthepeon.free.fr/ssemath/
    // for the exponential
    for (NeuronID n = 0 ; n < get_rank_size() ; ++n )
    {
		if (use_recovery)
		{
			handle_u[n] += dt/ms * izhi_a * ( u_coeffs[1] * handle_mem[n] + u_coeffs[0] - handle_u_temp[n] );
			handle_mem[n] += dt/ms/C*mV * ( mem_coeffs[2] * handle_mem[n]*handle_mem[n] + mem_coeffs[1] * handle_mem[n] + mem_coeffs[0] - handle_u_temp[n] + handle_inputcurrent[n]  ) ;
		}
		else
		{
			handle_mem[n] += dt/ms/C*mV * ( mem_coeffs[2] * handle_mem[n]*handle_mem[n] + mem_coeffs[1] * handle_mem[n] + mem_coeffs[0] - handle_u[n] + handle_inputcurrent[n]  ) ;
		}
    }

	// clip lower bound of membrane potential to account for bad (initial) conditions:
	//auryn_vector_float_clip( mem, V_min );
}

/// Integrate the internal state
/*!
       This method applies the Euler integration step to the membrane dynamics.
 */

void Izhikevich2003Group::check_peaks()
{
	for ( AurynState * i = mem->data ; i != mem->data+get_rank_size() ; ++i ) { // it's important to use rank_size here otherwise there might be spikes from units that do not exist
    	if ( *i > V_peak ) {
			NeuronID unit = i-mem->data;
			push_spike(unit);
		    set_val (mem, unit, V_reset); // reset
		    if (use_recovery)
		    	auryn_vector_float_add_constant(u,izhi_d); //refractory
	    		//auryn_vector_float_set(u,unit,d); //refractory
		}
	}

}

void Izhikevich2003Group::evolve()
{
	check_peaks(); // moved to front of function, so that the monitors can actually track the above-peak membrane potentials!

	double projMult = 10.0;
	auryn_vector_float_saxpy(projMult,g_ampa,t_exc);
	auryn_vector_float_set_zero(g_ampa);

	integrate_membrane_debug_two();

	auryn_vector_float_set_zero(t_exc);
}


AurynState Izhikevich2003Group::get_t_exc(NeuronID i)
{
	return get_val(t_exc,i);
}
AurynState Izhikevich2003Group::get_t_inh(NeuronID i)
{
	return get_val(t_inh,i);
}
AurynState Izhikevich2003Group::get_v_temp(NeuronID i)
{
	return get_val(v_temp,i);
}
AurynState Izhikevich2003Group::get_u_temp(NeuronID i)
{
	return get_val(u_temp,i);
}
AurynState Izhikevich2003Group::get_u(NeuronID i)
{
	return get_val(u,i);
}
AurynState Izhikevich2003Group::get_tempMemState(int i)
{
	return (AurynState)tempMemStates[i];
}



