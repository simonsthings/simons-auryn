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

#include <Izhikevich2007Group.h>


Izhikevich2007Group::Izhikevich2007Group( NeuronID size, AurynFloat load, NeuronID total ) : NeuronGroup(size,load,total)
{
	sys->register_spiking_group(this);
	if ( evolve_locally() ) init();
}


void Izhikevich2007Group::init()
{
	// BEGIN old stuff
	//e_rest = -70e-3;
	//e_reset = -70e-3;
	//set_ampa_nmda_ratio(1.0);
	t_leak = get_state_vector("t_leak");
	t_exc =  get_state_vector("t_exc");
	t_inh = get_state_vector("t_inh");
	// END old stuff

	
	// new stuff!
	if (true)
	{
		// To allow RS neuron of p. 274 of Izhikevich book, use these params:
		C = 100*pF;  // pF
		k = 0.7;                // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_thres = -40*mV;  // instantaneous threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_rest  = -60*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
		a = 0.03;
		b = -2;
		c = -50*mV; // mV
		d = 100;
		V_peak = 35*mV; //mV;
		// I_syn = 70; // pA
	}
	else
	{
		// To allow (0.04v^2 + 5v + 140) as in Izhikevich (2003) paper, use these params:
		C = 10*pF;  // pF
		k = 0.04;                // scaling factor k as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_thres = -42.3444*mV;  // instantaneous threshold potential v_t as in Izhikevich book p. 273 (Chapter 8.1.4)
		v_rest  = -82.6556*mV;   // resting potential v_r as in Izhikevich book p. 273 (Chapter 8.1.4)
		a = 0.02;
		b = 0.2;
		c = -40*mV; // mV
		d = 2;
		V_peak = 30*mV;
		// I_syn = 70; // pA
	}

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
	use_recovery = false;
	//e_rest = c;//*1e-3;  // mV. Not sure about scaling here, though.
	//e_reset = c;//*1e-3; // mV. Not sure about scaling here, though.
	doProperChannelStuff = false;   // for debugging.




	/** Begin set up Handles **/
    handle_mem = auryn_vector_float_ptr ( mem , 0 );
    handle_u = auryn_vector_float_ptr ( u , 0 );
    handle_u_temp = auryn_vector_float_ptr ( u_temp , 0 );
    handle_e_input = auryn_vector_float_ptr ( t_exc , 0 );
    /** End set up Handles **/

    // scale constants to keep time and other units SI-conform:
	scale_mem = dt/ms /(C/pF) * mV;
	scale_u = dt/ms;

	//cout << setiosflags(ios::scientific) << setprecision(6);
	cout << setiosflags(ios::fixed) << setprecision(6);

	cout << "scale_mem: " << scale_mem << endl;
	cout << "scale_u: " << scale_u << endl;

	cout << "v_thres: " << v_thres << endl;
	cout << "v_rest: " << v_rest << endl;

	clear();
}

void Izhikevich2007Group::clear()
{
	clear_spikes();
	for (NeuronID i = 0; i < get_rank_size(); ++i) {
	   //auryn_vector_float_set (mem, i, e_rest);
	   auryn_vector_float_set (mem, i, c);
	   auryn_vector_float_set (thr, i, 0.);
	   auryn_vector_float_set (g_ampa, i, 0.);
	   auryn_vector_float_set (g_gaba, i, 0.);
	   auryn_vector_float_set (g_nmda, i, 0.);

	   auryn_vector_float_set (u, i, 0.2*c);
	   auryn_vector_float_set (u_temp, i, 0.2*c);  // not really needed
	   auryn_vector_float_set (inputCurrents, i, 0.);
	   auryn_vector_float_set (backgroundCurrents, i, 0.);
	}
}

void Izhikevich2007Group::free() {
}

Izhikevich2007Group::~Izhikevich2007Group()
{
	if ( evolve_locally() ) free();
}


/// Integrate the internal state
void Izhikevich2007Group::integrate_membrane()
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
			handle_u[n] += scale_u * a * ( b * (handle_mem[n]-v_rest) - handle_u[n] );
			handle_mem[n] += scale_mem * ( k * (handle_mem[n]-v_rest) * (handle_mem[n]-v_thres) - handle_u_temp[n] + handle_inputcurrent[n] ) ;
		}
		else
		{
			cout << setiosflags(ios::fixed) << setprecision(6);
			cout << "mem: " << (handle_mem[n])*1e3 << endl;
			cout << "inputcurrent: " << (handle_inputcurrent[n])*1e3 << endl;
			cout << "diff mem-v_rest: " << (handle_mem[n]-v_rest)*1e3 << endl;
			cout << "diff mem-v_thresh: " << (handle_mem[n]-v_thres)*1e3 << endl << endl;

			// uses handle_u directly, instead of the copied handle_u_temp:
			handle_mem[n] += scale_mem * ( k * (handle_mem[n]-v_rest)/mV * (handle_mem[n]-v_thres)/mV - handle_u[n] + handle_inputcurrent[n] ) ;
			//handle_mem[n] += scale_mem * ( k * (handle_mem[n]-v_rest)/mV * (handle_mem[n]-v_thres)/mV - handle_u[n] + 0 ) ;
		}
    }


	// clip lower bound of membrane potential to account for bad (initial) conditions:
	//auryn_vector_float_clip( mem, V_min );
	tempMemStates[12] = handle_mem[0];


}

void Izhikevich2007Group::check_peaks()
{
	for ( AurynState * i = mem->data ; i != mem->data+get_rank_size() ; ++i ) { // it's important to use rank_size here otherwise there might be spikes from units that do not exist
    	if ( *i > V_peak ) {
			NeuronID unit = i-mem->data;
			push_spike(unit);
		    set_val (mem, unit, c); // reset
		    if (use_recovery)
		    	set_val (u, unit, d); //refractory
		}
	}

}

void Izhikevich2007Group::evolve()
{
	check_peaks(); // moved to front of function, so that the monitors can actually track the above-peak membrane potentials!


	//auryn_vector_float_copy(g_ampa,t_exc);
	//auryn_vector_float_set_zero(g_ampa);

	integrate_membrane();

	auryn_vector_float_set_zero(t_exc);

}



AurynState Izhikevich2007Group::get_t_exc(NeuronID i)
{
	return get_val(t_exc,i);
}

AurynState Izhikevich2007Group::get_t_inh(NeuronID i)
{
	return get_val(t_inh,i);
}

AurynState Izhikevich2007Group::get_v_temp(NeuronID i)
{
	return get_val(v_temp,i);
}

AurynState Izhikevich2007Group::get_u_temp(NeuronID i)
{
	return get_val(u_temp,i);
}

AurynState Izhikevich2007Group::get_u(NeuronID i)
{
	return get_val(u,i);
}

AurynState Izhikevich2007Group::get_tempMemState(int i)
{
	return (AurynState)tempMemStates[i];
}



