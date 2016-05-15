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
	use_recovery = true;

	select_example_parameter_set(Izhikevich_DefaultParametersets::default_2003);
	cout << "The Izhikevich coefficients are now: ";
	cout << " a2 = " << mem_coeffs[2];
	cout << ", a1 = " << mem_coeffs[1];
	cout << ", a0 = " << mem_coeffs[0] << endl;

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

void Izhikevich2003Group::select_example_parameter_set(Izhikevich_DefaultParametersets paramsetchoice)
{
	//AurynDouble k,v_rest,v_thres;
	V_min = -150; // -150 Volts. This is just a safeguard. TODO: Output warning if this is reached!

	switch (paramsetchoice)
	{
	case Izhikevich_DefaultParametersets::default_2003:
		a = 0.02;
		b = 0.2;
		c = -65; // mV, but don't use SI for now.
		d = 2;
		V_peak = 30;//*mV;
		// To allow (0.04v^2 + 5v + 140) as in Izhikevich (2003) paper, use these params:
		C = 1;
		k = 0.04;                   // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_thres = -42.3444;//*mV;  // soft threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_rest  = -82.6556;//*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
		//v_thres = 82.6556*mV;  // soft threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		//v_rest  = 42.3444*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
		// FIXME: not sure if the fit really is correct as v_rest and v_thres need to be positive for the coeffs to be computed correctly!
		V_reset = c;
		mem_coeffs[0] = 140;
		mem_coeffs[1] = 5;
		mem_coeffs[2] = 0.04 ;
		u_coeffs[0] = 0;
		u_coeffs[1] = b;
		break;
	case Izhikevich_DefaultParametersets::fitted_2003:
		a = 0.02;
		b = 0.2;
		c = -65;
		d = 2;
		V_peak = 30*mV;
		// To allow (0.04v^2 + 5v + 140) as in Izhikevich (2003) paper, use these params:
		C = 1;
		k = 0.04;                   // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
		//v_rest  = -82.6556;//*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
		//v_thres = -42.3444;//*mV;  // soft threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_rest  = 42.3444*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_thres = 82.6556*mV;  // soft threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		// FIXME: not sure if the fit really is correct as v_rest and v_thres need to be positive for the coeffs to be computed correctly!
		calculate_Izhikevich_coefficients(k,v_rest,v_thres,b);
		break;
	case Izhikevich_DefaultParametersets::A_2004:
		// To allow RS neuron of p. 274 of Izhikevich book, use these params:
		a = 0.03;
		b = -2;
		c = -50; // mV
		d = 100;
		V_peak = 35*mV;
		C = 100;  // pF
		k = 0.7;                // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_rest  = -60*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_thres = -40*mV;  // instantaneous threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		// I_syn = 70; // pA
		calculate_Izhikevich_coefficients(k,v_rest,v_thres,b);
		break;
	default:
		throw std::logic_error("Unhandled switch case!");
	}

}
void Izhikevich2003Group::calculate_Izhikevich_coefficients(AurynDouble k, AurynDouble v_rest, AurynDouble v_thres, AurynDouble b)
{
	mem_coeffs[0] = k * v_rest * v_thres;
	mem_coeffs[1] = k * (v_rest + v_thres);
	mem_coeffs[2] = k ;

	u_coeffs[0] = b*v_rest;
	u_coeffs[1] = b;

	V_reset = c*mV;
}
void Izhikevich2003Group::clear()
{
	clear_spikes();

	//auryn_vector_float_set (mem, i, e_rest);
	auryn_vector_float_set_all(mem, V_reset);
	auryn_vector_float_set_all (thr, 0.);
	auryn_vector_float_set_all (g_ampa, 0.);
	auryn_vector_float_set_all (g_gaba, 0.);
	auryn_vector_float_set_all (g_nmda, 0.);

	auryn_vector_float_set_all (u, b*c);
	auryn_vector_float_set_all (inputCurrents, 0.);
	auryn_vector_float_set_all (backgroundCurrents, 0.);
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
		// will refer to this later via handle_u_temp:
		auryn_vector_float_copy(u,u_temp);
	}

    // TODO we should vectorize this code and use some fast SSE
    // library such as http://gruntthepeon.free.fr/ssemath/
    // for the exponential
    for (NeuronID n = 0 ; n < get_rank_size() ; ++n )
    {
		if (use_recovery)
		{
			//double tempscaler = 1e3;
			double tempscaler = 1;
			cout << setiosflags(ios::fixed) << setprecision(6);
			cout << "with recovery!!" << endl;
			cout << "mem: " << (handle_mem[n])*tempscaler << endl;
			cout << "u:   " << (handle_u[n])  *tempscaler << endl;
			cout << "inputcurrent: " << (handle_inputcurrent[n])*tempscaler << endl;
			cout << "diff mem-v_rest: " << (handle_mem[n]-v_rest)*tempscaler << endl;
			cout << "diff mem-v_thresh: " << (handle_mem[n]-v_thres)*tempscaler << endl << endl;

			//handle_u[n] += dt/ms * a * ( b * (handle_mem[n]-v_rest)/mV - handle_u[n] );
			//handle_mem[n] += dt/ms * ( mem_coeffs[2] * handle_mem[n]*handle_mem[n] + mem_coeffs[1] * handle_mem[n] + mem_coeffs[0]  -  handle_u_temp[n] + handle_inputcurrent[n] ) ;

			handle_u[n] += dt * a * ( u_coeffs[1] * handle_mem[n] + u_coeffs[0] - handle_u_temp[n] );
			handle_mem[n] +=  dt * ( mem_coeffs[2] * handle_mem[n]*handle_mem[n] + mem_coeffs[1] * handle_mem[n] + mem_coeffs[0] - handle_u_temp[n] + handle_inputcurrent[n]  ) ;
		}
		else
		{

			//double tempscaler = 1e3;
			double tempscaler = 1;
			cout << setiosflags(ios::fixed) << setprecision(6);
			cout << "no recovery." << endl;
			cout << "mem: " << (handle_mem[n])*tempscaler << endl;
			cout << "u:   " << (handle_u[n])  *tempscaler << endl;
			cout << "inputcurrent: " << (handle_inputcurrent[n])*tempscaler << endl;
			cout << "diff mem-v_rest: " << (handle_mem[n]-v_rest)*tempscaler << endl;
			cout << "diff mem-v_thresh: " << (handle_mem[n]-v_thres)*tempscaler << endl << endl;


			// uses handle_u directly, instead of the copied handle_u_temp:
			//handle_mem[n] += scale_mem * ( k * (handle_mem[n]-v_rest)/mV * (handle_mem[n]-v_thres)/mV - handle_u[n] + 0 ) ;
			//handle_mem[n] += scale_mem * ( k * (handle_mem[n]-v_rest)/mV * (handle_mem[n]-v_thres)/mV - handle_u[n] + handle_inputcurrent[n] ) ;
		//handle_mem[n] += mV * dt * ( mem_coeffs[2] * handle_mem[n]*handle_mem[n]/mV/mV + mem_coeffs[1] * handle_mem[n]/mV + mem_coeffs[0]  -  -13 + handle_inputcurrent[n] ) ;
			//handle_mem[n] += dt/ms * ( mem_coeffs[2] * handle_mem[n]*handle_mem[n] + mem_coeffs[1] * handle_mem[n] + mem_coeffs[0]  -  handle_u[n] + 0 ) ;
			//handle_mem[n] += dt/ms * ( mem_coeffs[2] * handle_mem[n]*handle_mem[n] + mem_coeffs[1] * handle_mem[n] + mem_coeffs[0]  -  handle_u[n] + handle_inputcurrent[n] ) ;

			handle_mem[n] +=  dt * ( mem_coeffs[2] * handle_mem[n]*handle_mem[n] + mem_coeffs[1] * handle_mem[n] + mem_coeffs[0] - 13.0  + 000.0 + handle_inputcurrent[n]  ) ;
			//handle_mem[n] +=  dt * ( mem_coeffs[2] * handle_mem[n]*handle_mem[n] + mem_coeffs[1] * handle_mem[n] + mem_coeffs[0]  -  0 + 20 ) ;
		}
    }


	// clip lower bound of membrane potential to account for bad (initial) conditions:
	//auryn_vector_float_clip( mem, V_min );
	tempMemStates[12] = handle_mem[0];

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
	if (false)  // doProperChannelStuff
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

void Izhikevich2003Group::check_peaks()
{
	for ( AurynState * i = mem->data ; i != mem->data+get_rank_size() ; ++i ) { // it's important to use rank_size here otherwise there might be spikes from units that do not exist
    	if ( *i > V_peak ) {
			NeuronID unit = i-mem->data;
			push_spike(unit);
		    set_val (mem, unit, V_reset); // reset
		    if (use_recovery)
		    	set_val (u, unit, d); //refractory
		}
	}

}

void Izhikevich2003Group::evolve()
{
	check_peaks(); // moved to front of function, so that the monitors can actually track the above-peak membrane potentials!

	//auryn_vector_float_copy(g_ampa,t_exc);
	//auryn_vector_float_set_zero(g_ampa);

	//integrate_membrane();
	//integrate_membrane_debug();
	integrate_membrane_debug_two();

	// TODO: reactivate this?
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



