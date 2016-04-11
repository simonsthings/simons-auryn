/*
 * PolychronousPoissonGroup.cpp
 *
 *  Created on: 07.04.2016
 *      Author: simon
 */

#include "PolychronousPoissonGroup.h"

PolychronousPoissonGroup::PolychronousPoissonGroup(NeuronID N_total, NeuronID N_presenting, NeuronID N_subpresenting,
		AurynFloat duration, AurynFloat interval, NeuronID num_stimuli,  AurynDouble rate, string outputfilename) : PoissonGroup(N_total,rate)
{
	//PolychronousPoissonGroup::N_total = size;  // rename the size var for internal use
	init(N_presenting,N_subpresenting,duration,interval,num_stimuli,outputfilename);
}

PolychronousPoissonGroup::~PolychronousPoissonGroup()
{
	if ( evolve_locally() ) {
		patterntimesfile.close();
	}
}

void PolychronousPoissonGroup::init(NeuronID N_presenting, NeuronID N_subpresenting,
		AurynFloat duration, AurynFloat interval, NeuronID stimuli, string outputfilename)
{
	PolychronousPoissonGroup::N_presenting = N_presenting;
	PolychronousPoissonGroup::N_subpresenting = N_subpresenting;

	// Defaults:
	twoLegged = true;
	useRandomPermutations = false;
	participationProbability = 1.0;

	numPatterns = stimuli;
	patternDuration = duration/dt;
	patternInterval = interval/dt;
	logger->parameter("duration", (int)duration);
	logger->parameter("mean_isi", (int)patternInterval);
	max_patternDuration = patternDuration;
	max_patternInterval = patternInterval;

	representedValues.push_back(0.5);


	initBuffers(max_patternDuration);


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


	patterntimesfilename = outputfilename;
	if ( evolve_locally() && !outputfilename.empty()  ) {

		patterntimesfile.open(outputfilename.c_str(),ios::out);
		if (!patterntimesfile) {
			stringstream oss2;
			oss2 << "StructuredPoissonGroup:: Can't open output file " << outputfilename;
			logger->msg(oss2.str(),ERROR);
			exit(1);
		}
		patterntimesfile.setf(ios::fixed);
		patterntimesfile.precision(log(dt)/log(10)+1 );
	}






}

// Done:
void PolychronousPoissonGroup::initBuffers(int delaysteps)
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
		PolychronousPoissonGroup::constructNextBufferContents(theBuffer);
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
		PoissonGroup::evolve();
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
	auto myComparator = CompareIndicesByAnotherVectorValues<AurynTime>(unorderedLatencies);
	sort(orderingIndices->begin(), orderingIndices->end(), myComparator );
	//sort(orderingIndices->begin(), orderingIndices->end(), CompareIndicesByAnotherVectorValues<AurynTime>(unorderedLatencies) );
}


void PolychronousPoissonGroup::constructNextBufferContents(PermutableSpiketrainBuffer* theBuffer)
{
	int subpati = 0;

	//cout << "Constructing next buffer contents." << endl;

	return;

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
return;
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
		int countSU = 0;
		for (auto unit : spikingUnitIDs)
		{
			cout << "The spikingUnitID " << ++countSU << " is: " << unit << endl;
		}

		int countORD = 0;
		for (auto ordID : orderingIndices)
		{
			cout << "The orderingIndex " << ++countORD << " after sorting is: " << ordID << endl;
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

	theBuffer->permuteSpiketrains(&spikingUnitIDs,&orderingIndices);



//	theBuffer->clear();
//
//	for (int ordi = 0 ; ordi < orderingIndices.size() ; ++ordi)
//	{
//		theBuffer->setSingleSpike(spikingUnitIDs[ordi],unorderedLatencies[ordi]);
//	}



return;

	SpikeContainer newsubpatPresentingIDs;


	// permute buffer with these indices:




	// rearrange indices for twolegged or random permutation
	if (useRandomPermutations)
	{
		// rearrange permutation indices to represent the pattern for quantValue
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
	}




	// TODO: apply representedValues
	//unsigned int* permutationIndices;
	//theBuffer->permuteSpiketrains(permutationIndices);



}



void PolychronousPoissonGroup::evolve()
{
	//cout << "evolving PolychronousPoissonGroup!" << endl;

	if ( sys->get_clock() >= next_event ) {
		if ( stimulus_active ) {
			stimulus_active = false;
			//seed(sys->get_clock());
			next_event += (patternInterval-patternDuration);
			//cout << " ending stimulus! " << endl;
		} else {
			stimulus_active = true;
			current_stimulus = (current_stimulus+1)%numPatterns;
			x = 0;
			patterntimesfile << sys->get_time() << " " << current_stimulus << endl;
			//seed(current_stimulus+seedoffset);
			next_event += patternDuration;
			//cout << " starting new stimulus! " << endl;
		}
	}


	pregen_spikes.clear();

	if (stimulus_active)
	{
		if (!buffers[current_read_buffer].hasUnreadTimesteps())
		{
			// clear used buffer
			buffers[current_read_buffer].clear();


			// prepare used buffer for next use by filling with random data (possibly this could be done by separate thread in the future?)
			current_write_buffer = current_read_buffer;
			for (AurynTime ti = 0 ; ti < max_patternDuration ; ++ti)
			{
				// this could probably be sped up by direct implementation:
				pregen_spikes.clear();
				PoissonGroup::evolve();
				buffers[current_write_buffer].setSpikes(&pregen_spikes);
			}

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
		PoissonGroup::evolve();
		//pregen_spikes.clear(); // temporarily for debugging
	}

	// probably not the most efficient with the two+ for loops. But lets not optimise too early!
	for (NeuronID spike : pregen_spikes)
	{
		push_spike(spike);
	}

}







