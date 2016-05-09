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
	scale_ampa =  exp(-dt/tau_ampa) ;
	scale_gaba =  exp(-dt/tau_gaba) ;
	scale_thr = exp(-dt/tau_thr) ;
}

void Izhikevich2003Group::init()
{
	// BEGIN old stuff
	//e_rest = -70e-3;
	//e_reset = -70e-3;
	e_rev = -80e-3;
	thr_rest = -50e-3;
	dthr = 100e-3;
	tau_thr = 5e-3;
	tau_mem = 20e-3;
	tau_ampa = 1e-3;  // 5 ms
	tau_gaba = 10e-3;
	tau_nmda = 100e-3;
	set_ampa_nmda_ratio(1.0);
	calculate_scale_constants();
	t_leak = get_state_vector("t_leak");
	t_exc =  get_state_vector("t_exc");
	t_inh = get_state_vector("t_inh");
	// END old stuff
	
	// new stuff!
	mV = 1e-3;  // needed because Izhikevich neuron is originally defined in millivolts, not volts.
	ms = 1e-3;  // needed to allow for different integration time steps than 1 millisecond.
	a = 0.02;
	b = 0.2;
	c = -65;
	d = 2;
	//v = get_state_vector("v");  // using mem instead!
	u = get_state_vector("u");
	v_temp = get_state_vector("v_temp");
	v_temp2 = get_state_vector("v_temp2");
	u_temp = get_state_vector("u_temp");
	V_peak = 30*mV;
	V_min = -150; // -150 Volts. This is just a safeguard. TODO: Output warning if this is reached!
	inputCurrents =  get_state_vector("inputCurrents");
	backgroundCurrents = get_state_vector("backgroundCurrents");
	consistent_integration = true;
	use_recovery = true;
	e_rest = c*mV;//*1e-3;  // mV. Not sure about scaling here, though.
	e_reset = c*mV;//*1e-3; // mV. Not sure about scaling here, though.
	doProperChannelStuff = false;   // for debugging.

	// To allow (0.04v^2 + 5v + 140) as in Izhikevich (2003) paper, use these params:
	C = 1;
	k = 0.04;                   // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
	v_rest  = -82.6556;//*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
	v_thres = -42.3444;//*mV;  // soft threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)

	if (false)
	{
		// To allow RS neuron of p. 274 of Izhikevich book, use these params:
		C = 100;  // pF
		k = 0.7;                // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_rest  = -60;//*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_thres = -40;//*mV;  // instantaneous threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		a = 0.03;
		b = -2;
		c = -50; // mV
		d = 100;
		V_peak = 35; //mV;
		// I_syn = 70; // pA
	}

	clear();
}

void Izhikevich2003Group::clear()
{
	clear_spikes();
	for (NeuronID i = 0; i < get_rank_size(); ++i) {
	   //auryn_vector_float_set (mem, i, e_rest);
	   auryn_vector_float_set (mem, i, c*mV);
	   auryn_vector_float_set (thr, i, 0.);
	   auryn_vector_float_set (g_ampa, i, 0.);
	   auryn_vector_float_set (g_gaba, i, 0.);
	   auryn_vector_float_set (g_nmda, i, 0.);

	   auryn_vector_float_set (u, i, 0.2*e_rest);
	   auryn_vector_float_set (inputCurrents, i, 0.);
	   auryn_vector_float_set (backgroundCurrents, i, 0.);
	}
}

void Izhikevich2003Group::free() {
}

Izhikevich2003Group::~Izhikevich2003Group()
{
	if ( evolve_locally() ) free();
}

void Izhikevich2003Group::integrate_linear_nmda_synapses()
{
	// decay of ampa and gaba channel, i.e. multiply by exp(-dt/tau)
    auryn_vector_float_scale(scale_ampa,g_ampa);
    auryn_vector_float_scale(scale_gaba,g_gaba);

    // compute dg_nmda = (g_ampa-g_nmda)*dt/tau_nmda and add to g_nmda
	AurynFloat mul_nmda = dt/tau_nmda;
    auryn_vector_float_saxpy(mul_nmda,g_ampa,g_nmda);
	auryn_vector_float_saxpy(-mul_nmda,g_nmda,g_nmda);

    // excitatory
    auryn_vector_float_copy(g_ampa,t_exc);
    auryn_vector_float_scale(-A_ampa,t_exc);
    auryn_vector_float_saxpy(-A_nmda,g_nmda,t_exc);
    auryn_vector_float_mul(t_exc,mem);
    
    // inhibitory
    auryn_vector_float_copy(mem,t_inh);
    auryn_vector_float_add_constant(t_inh,-e_rev);
    auryn_vector_float_mul(t_inh,g_gaba);
}

/// Integrate the internal state
/*!
       This method applies the Euler integration step to the membrane dynamics.
 */
void Izhikevich2003Group::integrate_membrane_debug()
{
	//AurynFloat projMult = 2000;  // the projection multiplier as used in my matlab code. get rid of this at some point!

	AurynDouble h = dt/ms; //1;//0.1;
	AurynFloat V = mem->data[0];
	AurynFloat U = u->data[0]; // 0; //u->data[0];
	AurynFloat I_syn = t_exc->data[0];

	//AurynFloat I_scale = 0.02;  // before moving I_syn inside the brackets
	//AurynFloat I_scale = 0.005;   // still with dt = 0.1 ms, outside of bracket?

	AurynFloat I_scale;
	if (doProperChannelStuff)
	{
		if (dt == 1.0e-3)
			I_scale = 500; // dt = 1 ms,   Indirect. Works! (quite fine: some ghosts, but tuning to early spikes remains stable here! Ratemax=1000, rateZero=10s)
		else if (dt == 1.0e-4)
			//AurynFloat I_scale = 50;  // dt = 0.1 ms. Indirect. Works! (quite broad set of strong weights, but ok. Also, some ghost may be starting to materialise into a second peak.)
			I_scale = 50; // Indirect
		else
			cout << "Error error error error (doProperChannelStuff=true): dt=" << dt << endl;
	}
	else
	{
		if (dt == 1.0e-3)
			I_scale = 2;   // dt = 1 ms    Direct.   Works! (well, repeated learning, but its something! Ratemax=500, rateZero=11-12s)
		else if (dt == 1.0e-4)
			//AurynFloat I_scale = 20;  // dt = 0.1 ms. Direct.   Works! (stable red first 50 inputs, then darkblue for another 250 inputs. Then random-looking! :-o)
			I_scale = 20; // Direct
		else
			cout << "Error error error error (doProperChannelStuff=false): dt=" << dt << endl;
	}

	consistent_integration = true;

	// use standard forward Euler numerics in this case
	if (consistent_integration)
	{
		//mem->data[0] += mV * h * (0.04*V/mV*V/mV + 5.0*V/mV + 140 - U) + I_syn*I_scale;
		mem->data[0] += mV * h * (0.04*V/mV*V/mV + 5.0*V/mV + 140 - U + I_syn*I_scale) ;   // divided by C=1pF
	}
	else
	{
		mem->data[0] += 0.5*h * (0.04*V*V/mV + 5.0*V + 140*mV - U) + 0.5*I_syn*I_scale;
		AurynFloat V = mem->data[0];
		mem->data[0] += 0.5*h * (0.04*V*V/mV + 5.0*V + 140*mV - U) + 0.5*I_syn*I_scale;
	}

	if (use_recovery)
	{
		u->data[0] += h * a * (b * V - U);
	}

	tempMemStates[0] = h;
	tempMemStates[1] = V;
	tempMemStates[2] = U;
	tempMemStates[3] = I_syn;
	tempMemStates[4] = mem->data[0];

}
void Izhikevich2003Group::integrate_membrane()
{
	// use standard forward Euler numerics in this case
	if (consistent_integration)
	{
		/* BEGIN: Do the voltage */
		// the NEST code, for comparison:
		//S_.v_ += h*( 0.04*v_old*v_old + 5.0*v_old + 140.0 - u_old + S_.I_ + P_.I_e_)  +  B_.spikes_.get_value(lag) ;


		// need a temp vector:
		auryn_vector_float_copy(mem,v_temp);  // start with v_temp:=mem
		tempMemStates[0] = v_temp->data[0];

		// 0.04 * mem^2
		auryn_vector_float_mul(v_temp,mem); // mem ^ 2
		tempMemStates[1] = v_temp->data[0];
		auryn_vector_float_scale(0.04,v_temp); // 0.04 * [...]
		tempMemStates[2] = v_temp->data[0];

		auryn_vector_float_scale(1/mV,v_temp); // scale by 1mV, as done e.g. in the NeuroML example http://www.opensourcebrain.org/projects/izhikevichmodel/wiki/Wiki
		tempMemStates[3] = v_temp->data[0];


		// + 5 * mem
		auryn_vector_float_saxpy(5,mem,v_temp); // add 5*mem to v_temp
		tempMemStates[4] = v_temp->data[0];


		// start new temp variable. This is not actually needed, but increases readability for now.
		auryn_vector_float_set_all(v_temp2,140);  // + ( 140
		tempMemStates[5] = v_temp2->data[0];

		// add input currents. TODO: This line may need to be changed later!
		//auryn_vector_float_add(v_temp2,inputCurrents);

		// subtract u from 140:
		auryn_vector_float_saxpy(-1,u,v_temp2);  // - u
		tempMemStates[6] = v_temp2->data[0];

		// add constant background currents.
		//auryn_vector_float_add(v_temp,backgroundCurrents);

		// add any new spikes! Treat t_exc (?) as I_syn from neuroML model:
		auryn_vector_float_add(v_temp2,t_exc);  // + I_syn
		tempMemStates[7] = v_temp2->data[0];

		auryn_vector_float_scale(mV,v_temp2); // scale by 1mV, as done e.g. in the NeuroML example http://www.opensourcebrain.org/projects/izhikevichmodel/wiki/Wiki
		tempMemStates[8] = mem->data[0];

		// add the voltage-scaled last term to v_temp:
		auryn_vector_float_add(v_temp,v_temp2);  // + (140 - u + I_syn)
		tempMemStates[9] = v_temp->data[0];


		// multiply the whole thing with the simulation step size! TODO: check if I'm doing this right. This step may need to use "1" instead of dt. (In which case this line needs to be skipped!)
		auryn_vector_float_scale(dt/ms,v_temp);
		//auryn_vector_float_scale(ms/dt,v_temp); // likely wrong?
		tempMemStates[10] = v_temp->data[0];

		// add any new spikes!
		//auryn_vector_float_add(v_temp,t_exc);
		/* END: Do the voltage */



		/* BEGIN: Do the recovery */
		if (use_recovery)
		{
			// the NEST code, for comparison:
			//S_.u_ += h * P_.a_ * (P_.b_*v_old - u_old);

			// need a temp vector (bare with me: starting with contents of mem!):
			auryn_vector_float_copy(mem,u_temp); // fill temp vector from mem (=v)
			auryn_vector_float_scale(b,u_temp);  // b * v
			auryn_vector_float_saxpy(-1,u,u_temp);  // -1 * u + [...]

			auryn_vector_float_scale(a,u_temp);  // a * [...]
			// multiply the whole thing with the simulation step size! TODO: check if I'm doing this right. This step may need to use "1" instead of dt. (In which case this line needs to be skipped!)
			//auryn_vector_float_scale(dt,u_temp);
		}
		/* END: Do the recovery */

		/* BEGIN: Apply! */
		auryn_vector_float_add(mem,v_temp);
		tempMemStates[11] = mem->data[0];
		if (use_recovery)
			auryn_vector_float_add(u,u_temp);
		/* END: Apply! */
	}
	// use numerics published in Izhikevich (2003) in this case (not recommended)
	else
	{
//		S_.v_ += h/2.0 * ( 0.04*S_.v_*S_.v_ + 5.0*S_.v_ + 140.0 - S_.u_ + S_.I_ + P_.I_e_)
//	    		   +  B_.spikes_.get_value(lag);
//		S_.v_ += h/2.0 * ( 0.04*S_.v_*S_.v_ + 5.0*S_.v_ + 140.0 - S_.u_ + S_.I_ + P_.I_e_)
//	    		   +  B_.spikes_.get_value(lag);
//		S_.u_ += h * P_.a_*(P_.b_*S_.v_ - S_.u_);
	}

	// lower bound of membrane potential
	auryn_vector_float_clip( mem, V_min );
	tempMemStates[12] = mem->data[0];
}

//void IzhikevichGroup::integrate_membrane_IFNeuron()
//{
//	// moving threshold
//    auryn_vector_float_scale(scale_thr,thr);
//
//    // leak
//	auryn_vector_float_copy(mem,t_leak);
//    auryn_vector_float_add_constant(t_leak,-e_rest);
//
//    // membrane dynamics
//	AurynFloat mul_tau_mem = dt/tau_mem;
//    auryn_vector_float_saxpy(mul_tau_mem,t_exc,mem);
//    auryn_vector_float_saxpy(-mul_tau_mem,t_inh,mem);
//    auryn_vector_float_saxpy(-mul_tau_mem,t_leak,mem);
//}

void Izhikevich2003Group::check_peaks()
{
	for ( AurynState * i = mem->data ; i != mem->data+get_rank_size() ; ++i ) { // it's important to use rank_size here otherwise there might be spikes from units that do not exist
    	if ( *i > V_peak ) {
			NeuronID unit = i-mem->data;
			push_spike(unit);
		    set_val (mem, unit, e_reset); // reset
		    if (use_recovery)
		    	set_val (u, unit, d); //refractory
		}
	}

}

void Izhikevich2003Group::evolve()
{
	check_peaks(); // moved to front of function, so that the monitors can actually track the above-peak membrane potentials!

	if (doProperChannelStuff)
	{
		integrate_linear_nmda_synapses();
	}
	else
	{
		auryn_vector_float_copy(g_ampa,t_exc);
		auryn_vector_float_set_zero(g_ampa);
	}

	//integrate_membrane();
	integrate_membrane_debug();
}


void Izhikevich2003Group::set_tau_mem(AurynFloat taum)
{
	tau_mem = taum;
	calculate_scale_constants();
}

AurynFloat Izhikevich2003Group::get_tau_mem()
{
	return tau_mem;
}

void Izhikevich2003Group::set_tau_ampa(AurynFloat taum)
{
	tau_ampa = taum;
	calculate_scale_constants();
}

AurynFloat Izhikevich2003Group::get_tau_ampa()
{
	return tau_ampa;
}

void Izhikevich2003Group::set_tau_gaba(AurynFloat taum)
{
	tau_gaba = taum;
	calculate_scale_constants();
}

AurynFloat Izhikevich2003Group::get_tau_gaba()
{
	return tau_gaba;
}

void Izhikevich2003Group::set_tau_nmda(AurynFloat taum)
{
	tau_nmda = taum;
	calculate_scale_constants();
}

AurynFloat Izhikevich2003Group::get_tau_nmda()
{
	return tau_nmda;
}

void Izhikevich2003Group::set_ampa_nmda_ratio(AurynFloat ratio)
{
 	A_ampa = ratio/(ratio+1.0);
	A_nmda = 1./(ratio+1.0);
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



