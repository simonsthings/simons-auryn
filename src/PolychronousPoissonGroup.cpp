/*
 * PolychronousPoissonGroup.cpp
 *
 *  Created on: 07.04.2016
 *      Author: simon
 */

#include "PolychronousPoissonGroup.h"

using namespace auryn;

boost::mt19937 PolychronousPoissonGroup::gen_PPG = boost::mt19937();

PolychronousPoissonGroup::PolychronousPoissonGroup(NeuronID N_total, NeuronID N_presenting, NeuronID N_subpresenting,
												   AurynFloat duration, AurynFloat interval, NeuronID num_stimuli,  AurynDouble rate, string outputfilename) : SpikingGroup(N_total)
{
	//PolychronousPoissonGroup::N_total = size;  // rename the size var for internal use
	init(rate,N_presenting,N_subpresenting,duration,interval,num_stimuli,outputfilename);
}

PolychronousPoissonGroup::~PolychronousPoissonGroup()
{
	if ( evolve_locally() ) {
		patterntimesfile.close();
	}
}

void PolychronousPoissonGroup::init(AurynDouble rate, NeuronID N_presenting, NeuronID N_subpresenting,
									AurynFloat duration, AurynFloat interval, NeuronID stimuli, string outputfilename)
{
	auryn::sys->register_spiking_group(this);

	PolychronousPoissonGroup::N_presenting = N_presenting;
	PolychronousPoissonGroup::N_subpresenting = N_subpresenting;

	// Defaults:
	twoLegged = false;
	useRandomPermutations = false;
	participationProbability = 1.0;

	numPatterns = stimuli;
	patternDuration = (AurynTime) (duration/auryn_timestep);
	patternInterval = (AurynTime) (interval/auryn_timestep);
	logger->parameter("duration", (int)duration);
	logger->parameter("mean_isi", (int)patternInterval);
	//max_patternDuration = patternDuration;
	//max_patternInterval = patternInterval;

	representedValues.push_back(0.5);

	// state information:
	stimulus_active = false;
	current_stimulus = 0;
	next_event = 0;


	stringstream oss;
	oss << "PolychronousPoissonGroup:: Set up with stimulus_duration="
		<< patternDuration
		<< " and mean_isi="
		<< patternInterval;
	logger->msg(oss.str(),NOTIFICATION);

//	if ( evolve_locally() )
//	{

	logger->msg("Now setting up the random stuff...",NOTIFICATION);
	logger->parameter("numPatterns", (int)numPatterns);

	dist_PPG = new boost::uniform_01<>();

	die_PPG = new boost::variate_generator<boost::mt19937 &, boost::uniform_01<> >(PolychronousPoissonGroup::gen_PPG, *dist_PPG);
	salt_PPG = sys->get_seed();
	seed(sys->get_seed());
	x = 0;
	set_rate( rate );


	logger->msg("Finished setting up the random stuff.",NOTIFICATION);
//	}

	patterntimesfilename = outputfilename;
	if ( evolve_locally() && !outputfilename.empty()  ) {

		patterntimesfile.open(outputfilename.c_str(),ios::out);
		if (!patterntimesfile) {
			stringstream oss2;
			oss2 << "PolychronousPoissonGroup:: Can't open output file " << outputfilename;
			logger->msg(oss2.str(),ERROR);
			exit(1);
		}
		patterntimesfile.setf(ios::fixed);
		patterntimesfile.precision(log(auryn_timestep)/log(10)+1 );
	}

	// should be at the end of the init function:
	initBuffers(patternDuration);
}

void PolychronousPoissonGroup::seed(unsigned int s)
{
	std::stringstream oss;
	oss << "PolychronousPoissonGroup:: Seeding rank " << sys->mpi_rank() << " with seed " << s << " and salt " << salt_PPG;
	logger->msg(oss.str(),NOTIFICATION);

	gen_PPG.seed( s + salt_PPG );
}



// Done:
void PolychronousPoissonGroup::initBuffers(AurynTime delaysteps)
{
	buffers.clear();

	int numBuffers = 3;
	// allocate memory and generate noise:
	for (int bi = 0; bi < numBuffers; ++bi)
	{
		PermutableSpiketrainBuffer* newbuffer = new PermutableSpiketrainBuffer(size,delaysteps);
		generateNoiseAhead(newbuffer, delaysteps);
		buffers.push_back(*newbuffer);
	}
	// fill with patterns:
	for (int bi = 0; bi < (numBuffers); ++bi)
	{
		PermutableSpiketrainBuffer* theBuffer = &buffers[bi];
		constructNextBufferContents(theBuffer);
	}
	current_read_buffer = 1;
	current_write_buffer = 1;
}

// Done:
void PolychronousPoissonGroup::generateNoiseAhead(PermutableSpiketrainBuffer* buffer, AurynTime steps)
{
	buffer->clear();
	//for (PermutableSpiketrainBuffer spikes_at_one_timestep : (*buffer))
	for (AurynTime ti = 0 ; ti < steps ; ++ti)
	{
		pregen_spikes.clear();
		PGevolve();
		buffer->setSpikes(&pregen_spikes);
	}
}

/**
 * Overwrite the location to which the freshly generated spike will be sent!
 * This allows us to use the evolve function of PoissonGroup to pre-generate
 * many timesteps worth of spikes, and also to temporarily store them differently
 * than by using an array of SpikeContainers (as is done in the SpikeDelay class).
 * @param theSpikingNeuron
 */
void PolychronousPoissonGroup::distribute_spike(NeuronID theSpikingNeuron)
{
	pregen_spikes.push_back ( theSpikingNeuron );
}

void PolychronousPoissonGroup::argsort(vector<AurynTime>* unorderedLatencies, SpikeContainer* orderingIndices)
{
	// initialize original index locations
	orderingIndices->resize(unorderedLatencies->size());
	for (size_t i = 0; i != orderingIndices->size(); ++i) (*orderingIndices)[i] = i;
	//std::iota( orderingIndices->begin(), orderingIndices->end(), 0 );


//	int countORD = 0;
//	for (auto ordID : (*orderingIndices))
//	{
//	cout << "In argsort(), the orderingIndex " << ++countORD << " is: " << ordID << endl;
//	}


	//const vector<AurynTime> &v = (*unorderedLatencies);
	// sort indexes based on comparing values in v
	//sort(orderingIndices->begin(), orderingIndices->end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
	const CompareIndicesByAnotherVectorValues <AurynTime> &myComparator = CompareIndicesByAnotherVectorValues<AurynTime>(unorderedLatencies);
	std::sort(orderingIndices->begin(), orderingIndices->end(), myComparator );
	//sort(orderingIndices->begin(), orderingIndices->end(), CompareIndicesByAnotherVectorValues<AurynTime>(unorderedLatencies) );
}


void PolychronousPoissonGroup::constructNextBufferContents(PermutableSpiketrainBuffer* theBuffer)
{
	int subpati = 0;

	//cout << "Constructing next buffer contents." << endl;

	// get the quantized value for this varchoices value. Round it so that we have recurring patterns:
	int quantValue = int(round(N_presenting * representedValues[subpati] ));

	spikingUnitIDs.clear();
	theBuffer->getUnitsThatSpikeAtLeastOnce(&spikingUnitIDs,N_subpresenting);

	// apply participationProbability:
    if (participationProbability < 1.0)
    {
    	for (auto iter = spikingUnitIDs.begin() ; iter != spikingUnitIDs.end() ; ++iter )
    	{
    	    if (participationProbability < float_dist(poly_gen))
    	    	spikingUnitIDs.erase(iter);
    	}
    }

//	int countSU = 0;
//	for (auto unit : spikingUnitIDs)
//	{
//		cout << "The spikingUnitID " << ++countSU << " is: " << unit << endl;
//	}

    // get some latencies to sort on:
    unorderedLatencies.clear();
	theBuffer->getRandomSpiketimeOfActiveUnits(&unorderedLatencies,&spikingUnitIDs);
	//theBuffer->getActiveUnitsAndRandomSpiketimes(&spikingUnitIDs,&unorderedLatencies,N_subpresenting,poly_gen,int_dist);
	//theBuffer->getRandomSpiketimeOfActiveUnits(&unorderedLatencies, &spikingUnitIDs, poly_gen, int_dist);

	//cout << "PolychronousPoissonGroup::constructNextBufferContents has been executed as far as implemented." << endl;

	// get sorting indices:
//	int count = 0;
//	for (auto lat : unorderedLatencies)
//	{
//		cout << "The lat " << ++count << " is: " << lat << endl;
//	}

	vector<NeuronID> orderingIndices(spikingUnitIDs.size());

	argsort(&unorderedLatencies, &orderingIndices);

	if (false)
	{
		//int countSU = 0;
		for (int countSU = 0 ; countSU < spikingUnitIDs.size() ; ++countSU)
		//for (auto unit : spikingUnitIDs)
		{
			cout << "The spikingUnitID " << countSU << " is: " << spikingUnitIDs[countSU] << endl;
		}

		//int countORD = 0;
		for (int countORD = 0 ; countORD < orderingIndices.size() ; ++countORD)
		//for (auto ordID : orderingIndices)
		{
			cout << "The orderingIndex " << countORD << " after sorting is: " << orderingIndices[countORD] << endl;
		}


		for (int ordi = 0 ; ordi < orderingIndices.size() ; ++ordi)
		{
			NeuronID theOldLatency = unorderedLatencies[ ordi ];
			cout << "The latency " << ordi << " before sorting is: " << theOldLatency << endl;
		}

		for (int ordi = 0 ; ordi < orderingIndices.size() ; ++ordi)
		{
			NeuronID theNewLatency = unorderedLatencies[ orderingIndices[ordi] ];
			cout << "The latency " << ordi << " after sorting is: " << theNewLatency << endl;
		}
	}

	// TODO: remove this after implementing below:
	theBuffer->permuteSpiketrains(&spikingUnitIDs,&orderingIndices);



//	theBuffer->clear();
//
//	for (int ordi = 0 ; ordi < orderingIndices.size() ; ++ordi)
//	{
//		theBuffer->setSingleSpike(spikingUnitIDs[ordi],unorderedLatencies[ordi]);
//	}



return;

	// make array of unit indices and fill with ID. (see later if this can be optimised)
	vector<NeuronID> newsubpatPresentingIDs(N_presenting);
	for (size_t i = 0; i != newsubpatPresentingIDs.size(); ++i) (newsubpatPresentingIDs)[i] = i;
	//std::iota( orderingIndices->begin(), orderingIndices->end(), 0 );


	// rearrange indices for twolegged or random permutation and apply value:
	if (useRandomPermutations)
	{
		// rearrange permutation indices to represent the pattern for quantValue
		// permute buffer with these indices:
		//vector<NeuronID> permutationIndices(spikingUnitIDs.size());
	}
	else
	{
		if (twoLegged)
		{
			// make two legs
		}
		else
		{
			// leave as-is
		}

		// now shift by value.
	}




	// TODO: apply representedValues
	//unsigned int* permutationIndices;
	//theBuffer->permuteSpiketrains(permutationIndices);

	theBuffer->permuteSpiketrains(&spikingUnitIDs,&orderingIndices);


}



void PolychronousPoissonGroup::evolve()
{
	//cout << "evolving PolychronousPoissonGroup!" << endl;

	// each stimulus only starts for one clock cylce. At all other times, the size is zero.
	stimuli_immediate.clear();

	if ( sys->get_clock() >= next_event ) {
		if ( stimulus_active ) {
			// switch off stimulus:
			stimulus_active = false;
			//seed(sys->get_clock());
			next_event += (patternInterval-patternDuration);
			//cout << " ending stimulus! " << endl;
		} else {
			// switch on stimulus:
			stimulus_active = true;
			current_stimulus = (current_stimulus+1)%numPatterns;
			stimuli_immediate.push_back(current_stimulus);
			x = 0;
			patterntimesfile << sys->get_time() << " " << current_stimulus << endl;
			//seed(current_stimulus+seedoffset);
			next_pattern_onset = next_event + patternInterval;
			next_event += patternDuration;
			//cout << " starting new stimulus! " << endl;
		}
	}


	pregen_spikes.clear();

	if (stimulus_active)
	{
		if (!buffers[current_read_buffer].hasUnreadTimesteps())
		{
			// using generateNoiseAhead() now instead of this:
			// clear used buffer
			//buffers[current_read_buffer].clear();

			// the last read buffer is now our next write buffer ;)
			current_write_buffer = current_read_buffer;

			// prepare used buffer for next use by filling with random data (possibly this could be done by separate thread in the future?)
			generateNoiseAhead(&buffers[current_write_buffer], patternDuration);

			// using generateNoiseAhead() now instead of this:
//			for (AurynTime ti = 0 ; ti < max_patternDuration ; ++ti)
//			{
//				// this could probably be sped up by direct implementation:
//				pregen_spikes.clear();
//				PGevolve();
//				buffers[current_write_buffer].setSpikes(&pregen_spikes);
//			}

			// then switch buffers
			current_read_buffer = (current_read_buffer + 1)%buffers.size();

			// generate pattern in new buffer
			PolychronousPoissonGroup::constructNextBufferContents(&buffers[current_read_buffer]);
		}
		// take from active buffer
		buffers[current_read_buffer].getSpikes(&pregen_spikes);
	}
	else
	{
		// just generate newrandom data:
		PGevolve();
		//pregen_spikes.clear(); // temporarily for debugging
	}

	// probably not the most efficient with the two+ for loops. But lets not optimise too early!
	for (NeuronID spike : pregen_spikes)
	{
		push_spike(spike);
	}


	checkAndUpdateTestingProtocol();
}


void PolychronousPoissonGroup::set_rate(AurynDouble  rate)
{
	cout << "The rate is being set to " << rate << "Hz." << endl;
	lambda_PPG = 1.0/(1.0/rate-auryn_timestep);
	if ( evolve_locally() ) {
		if ( rate > 0.0 ) {
			AurynDouble r = -log((*die_PPG)()+1e-128)/lambda_PPG;
			x = (NeuronID)(r/auryn_timestep+0.5);
		} else {
			// if the rate is zero this triggers one spike at the end of time/groupsize
			// this is the easiest way to take care of the zero rate case, which should
			// be avoided in any case.
			x = std::numeric_limits<NeuronID>::max();
		}
	}
	cout << "In accordance with the sim time step, lambda_PPG has become " << lambda_PPG << " ." << endl;
}
AurynDouble  PolychronousPoissonGroup::get_rate()
{
	return lambda_PPG;
}

vector<PatternID> PolychronousPoissonGroup::get_stimuli_immediate()
{
	return stimuli_immediate;
}

PatternID PolychronousPoissonGroup::getNumPatterns()
{
	return numPatterns;
}


void PolychronousPoissonGroup::setTestingProtocol(vector<AurynFloat> theTestprotocolDurations, vector<AurynFloat> theTestprotocolPatternintervals)
{
	testprotocolDurations = theTestprotocolDurations;
	testprotocolPatternintervals = theTestprotocolPatternintervals;

	// assumes that testprotocolDurations contains at least one element. Todo: ensure testprotocolDurations always contains at least one element!
	current_phase_id = 0;
	patternInterval = (AurynTime)(testprotocolPatternintervals[current_phase_id] / auryn_timestep);
	next_phase_start_clock = (AurynTime)(testprotocolDurations[current_phase_id] / auryn_timestep);
}

/**
 * Shift through the predefined phases of the testing protocol by adjusting the patternInterval accordingly.
 *
 * Note: As long as the total simulation time is the summation of all testprotocolDurations, this should not run into problems:
 * next_phase_id will then always be smaller than the length of testprotocolDurations.
 */
void PolychronousPoissonGroup::checkAndUpdateTestingProtocol()
{

	if ( sys->get_clock() > next_phase_start_clock )
	{
		current_phase_id++;
		patternInterval = (AurynTime)(testprotocolPatternintervals[current_phase_id] / auryn_timestep);
		next_phase_start_clock += (AurynTime)(testprotocolDurations[current_phase_id] / auryn_timestep);
	}

}




void PolychronousPoissonGroup::PGevolve()
{
	while ( x < get_rank_size() ) {
		//cout << " n" << x;
		distribute_spike ( x );

		AurynDouble r = -log((*die_PPG)()+1e-128)/lambda_PPG;
		//cout << " r" << r;
		//cout << " l" << lambda_PPG;
//		AurynDouble r = 20;
		// we add 1.5: one to avoid two spikes per bin and 0.5 to
		// compensate for rounding effects from casting
		x += (NeuronID)(r/auryn_timestep+1.5);
		// beware one induces systematic error that becomes substantial at high rates, but keeps neuron from spiking twice per time-step
	}
	x -= get_rank_size();
}

AurynTime PolychronousPoissonGroup::get_next_pattern_onset() const
{
	return next_pattern_onset;
}
































