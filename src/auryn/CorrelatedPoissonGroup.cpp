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
*
* If you are using Auryn or parts of it for your work please cite:
* Zenke, F. and Gerstner, W., 2014. Limits to high-speed simulations 
* of spiking neural networks using general-purpose computers. 
* Front Neuroinform 8, 76. doi: 10.3389/fninf.2014.00076
*/

#include "CorrelatedPoissonGroup.h"

using namespace auryn;

boost::mt19937 CorrelatedPoissonGroup::gen = boost::mt19937(); 

void CorrelatedPoissonGroup::init(AurynDouble  rate, NeuronID gsize, AurynDouble timedelay )
{
	auryn::sys->register_spiking_group(this);
	if ( evolve_locally() ) {
		lambda = rate;

		groupsize = global2rank(gsize);
		ngroups = size/gsize;
		remainersize = get_rank_size()-ngroups*groupsize;
		delay = timedelay/auryn_timestep;
		offset = 0;

		amplitude = 1.0;
		target_amplitude = amplitude;
		tau_amplitude = 60;
		thr = 1e-6;
		mean = 0.0;

		timescale = 50e-3;
		set_stoptime(0);

		std::stringstream oss;
		oss << "CorrelatedPoissonGroup:: Initializing with " 
			<< ngroups << " of size " 
			<< groupsize << " ( "
			<< " delay=" << timedelay << ", "
			<< " amplitude=" << amplitude << ", "
			<< " timescale=" << timescale 
			<< " )";
		auryn::logger->msg(oss.str(),NOTIFICATION);

		// delay_o = new AurynDouble [delay*(ngroups-1)];
		o = 1.0;
		delay_o = new AurynDouble [delay*ngroups];
		for ( unsigned int i = 0 ; i < delay*ngroups ; ++i )
			delay_o[i] = 1.0;

		dist = new boost::uniform_01<> ();
		die  = new boost::variate_generator<boost::mt19937&, boost::uniform_01<> > ( gen, *dist );
		seed(sys->mpi_rank()); // seeding problem

		x = new NeuronID [ngroups];
		for ( unsigned int i = 0 ; i < ngroups ; ++i ) {
			AurynDouble r = log(1-(AurynDouble)(*die)())/(-lambda);
			x[i] += (NeuronID)(r/auryn_timestep); 
		}

		oss.str("");
		oss << "CorrelatedPoissonGroup:: Seeding with " << sys->mpi_rank();
		auryn::logger->msg(oss.str(),NOTIFICATION);
	}
}

CorrelatedPoissonGroup::CorrelatedPoissonGroup(NeuronID n, 
		AurynDouble  rate, 
		NeuronID gsize, 
		AurynDouble timedelay ) : SpikingGroup( n , CORRELATEDPOISSON_LOAD_MULTIPLIER*rate ) 
{
	init(rate, gsize,  timedelay);
}

CorrelatedPoissonGroup::~CorrelatedPoissonGroup()
{
	if ( evolve_locally() ) {
		delete dist;
		delete die;
		delete delay_o;
	}
}

void CorrelatedPoissonGroup::set_rate(AurynDouble  rate)
{
	lambda = rate;
}

void CorrelatedPoissonGroup::set_threshold(AurynDouble  threshold)
{
	thr = std::max(1e-6,threshold);
}

AurynDouble  CorrelatedPoissonGroup::get_rate()
{
	return lambda;
}


void CorrelatedPoissonGroup::evolve()
{
	// check if the group has timed out
	if ( tstop && auryn::sys->get_clock() > tstop ) return;

	// move amplitude
	amplitude += (target_amplitude-amplitude)*auryn_timestep/tau_amplitude;

	// integrate Ornstein Uhlenbeck process
	o += ( mean - o )*auryn_timestep/timescale;

	// noise increment
	o += 2.0*((AurynDouble)(*die)()-0.5)*sqrt(auryn_timestep/timescale)*amplitude;

	int len = delay*ngroups;
	delay_o[auryn::sys->get_clock()%len] = std::max(thr,o*lambda);

	for ( unsigned int g = 0 ; g < ngroups ; ++g ) {
		AurynDouble grouprate = delay_o[(auryn::sys->get_clock()-(g+offset)*delay)%len];
		AurynDouble r = -log(1-(AurynDouble)(*die)())/(auryn_timestep*grouprate); // think before tempering with this! 
		// I already broke the corde here once!
		x[g] = (NeuronID)(r); 
		while ( x[g] < groupsize ) {
			push_spike ( g*groupsize + x[g] );
			AurynDouble r = -log(1-(AurynDouble)(*die)())/(auryn_timestep*grouprate);
			x[g] += (NeuronID)(r); 
		}
	}
}

void CorrelatedPoissonGroup::seed(int s)
{
		gen.seed(s); 
}

void CorrelatedPoissonGroup::set_amplitude(AurynDouble amp)
{
	amplitude = amp;
	auryn::logger->parameter("amplitude",amplitude);
}

void CorrelatedPoissonGroup::set_target_amplitude(AurynDouble amp)
{
	target_amplitude = amp;
	auryn::logger->parameter("target_amplitude",target_amplitude);
}

void CorrelatedPoissonGroup::set_timescale(AurynDouble scale)
{
	timescale = scale;
	auryn::logger->parameter("timescale",timescale);
}

void CorrelatedPoissonGroup::set_tau_amplitude(AurynDouble tau)
{
	tau_amplitude = tau;
	auryn::logger->parameter("tau_amplitude",tau_amplitude);
}

void CorrelatedPoissonGroup::set_offset(int off)
{
	offset = off;
	auryn::logger->parameter("offset",offset);
}

void CorrelatedPoissonGroup::set_stoptime(AurynDouble stoptime)
{
	tstop = stoptime*auryn_timestep;
	auryn::logger->parameter("stoptime",stoptime);
}
