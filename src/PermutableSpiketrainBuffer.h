/*
 * PermutableSpiketrainBufferNew.h
 *
 *  Created on: 07.04.2016
 *      Author: simon
 */

#ifndef SRC_PERMUTABLESPIKETRAINBUFFER_H_
#define SRC_PERMUTABLESPIKETRAINBUFFER_H_

#include "auryn_definitions.h"
#include "System.h"
#include "SpikingGroup.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <random>

using namespace std;

/**
 * A quick implementation to easily permutate spike trains. This likely becomes inefficient for
 * very large numbers of neurons, but should be fine for 10000 neurons or less.
 */
class PermutableSpiketrainBuffer
{
protected:
	vector< bool* > theSpiketrains; ///< lets see if nesting vectors is efficient enough. Don't case for now.
	NeuronID N_total;
	AurynTime timesteps_total;
	bool ispermuted;
	SpikeContainer permutationIndices;

	unsigned int write_position;
	unsigned int read_position;
	vector<AurynTime> temp_spikepositions;

    boost::random::mt19937 temp_gen;
    boost::random::uniform_int_distribution<> int_dist;


//public:
protected:
	void setSingleSpike(NeuronID nid, AurynTime timestep);
	void clearSingleSpike(NeuronID nid, AurynTime timestep);
	bool hasSpikeAt(NeuronID nid, AurynTime timestep);

	void setSpikes(SpikeContainer* spikes, AurynTime timestep);
	void getSpikesAtTimestep(SpikeContainer* targetspikecontainer, AurynTime timestep);

public:
	void setSpikes(SpikeContainer* spikes);
	void getSpikes(SpikeContainer* targetspikecontainer);
	bool hasUnreadTimesteps();
	void clear();
	void permuteSpiketrains(SpikeContainer* whoToPermute, SpikeContainer* thePermutationIndices);
	void getUnitsThatSpikeAtLeastOnce(SpikeContainer* targetspikecontainer, NeuronID N_presenting);
	void getRandomSpiketimeOfActiveUnits(vector<AurynTime>* randomspiketimes, const SpikeContainer* activeUnits);
	//void getActiveUnitsAndRandomSpiketimes(SpikeContainer* activeUnits, vector<AurynTime>* randomspiketimes,  NeuronID N_presenting,  boost::random::mt19937 poly_gen,  boost::random::uniform_int_distribution<> int_dist);


	PermutableSpiketrainBuffer(NeuronID N_total, AurynTime timesteps_total);
	virtual ~PermutableSpiketrainBuffer();
};

#endif /* SRC_PERMUTABLESPIKETRAINBUFFER_H_ */
