/*
 * PermutableSpiketrainBufferNew.cpp
 *
 *  Created on: 07.04.2016
 *      Author: simon
 */

#include <PermutableSpiketrainBuffer.h>


PermutableSpiketrainBuffer::PermutableSpiketrainBuffer(NeuronID N_total, AurynTime timesteps_total)
{
	ispermuted = false;
	PermutableSpiketrainBuffer::N_total = N_total;
	PermutableSpiketrainBuffer::timesteps_total = timesteps_total;
	write_position = 0;
	read_position = 0;

	//permutationIndices.reserve(N_total);
	permutationIndices.resize(N_total);
	theSpiketrains.reserve(N_total);
	for (int ni = 0 ; ni < N_total ; ++ni)
	{
		bool* spiketrain = new bool[timesteps_total];
		theSpiketrains.push_back(spiketrain);
	}
	clear();

	//int_dist = new boost::random::uniform_int_distribution<>(0,100);
}
PermutableSpiketrainBuffer::~PermutableSpiketrainBuffer()
{
	// deletion of theSpiketrains contents done by vector class.
}


void PermutableSpiketrainBuffer::setSpikes(SpikeContainer* spikes, AurynTime timestep)
{
	for (NeuronID nid : (*spikes))
	{
		theSpiketrains[nid][timestep] = true;
	}
}

void PermutableSpiketrainBuffer::getSpikesAtTimestep(SpikeContainer* targetspikecontainer, AurynTime timestep)
{
	for (NeuronID ni = 0 ; ni < N_total ; ++ni)
	{
		auto permutedNeuron = permutationIndices[ni];
		if (theSpiketrains[permutedNeuron][timestep])
		{
			targetspikecontainer->push_back(permutedNeuron);
		}
	}
}






void PermutableSpiketrainBuffer::setSingleSpike(NeuronID nid, AurynTime timestep)
{
	theSpiketrains[nid][timestep] = true;
}

void PermutableSpiketrainBuffer::clearSingleSpike(NeuronID nid, AurynTime timestep)
{
	theSpiketrains[nid][timestep] = false;
}

bool PermutableSpiketrainBuffer::hasSpikeAt(NeuronID nid, AurynTime timestep)
{
	return theSpiketrains[nid][timestep];
}





void PermutableSpiketrainBuffer::setSpikes(SpikeContainer* spikes)
{
	if (write_position < timesteps_total)
		setSpikes(spikes,write_position++);
	// else do nothing.
}
void PermutableSpiketrainBuffer::getSpikes(SpikeContainer* targetspikecontainer)
{
	if (read_position < write_position)
		getSpikesAtTimestep(targetspikecontainer,read_position++);
	// else do not fill anything.
}
bool PermutableSpiketrainBuffer::hasUnreadTimesteps()
{
	return (read_position < write_position);
}
/**
 * TODO: Make a lot more efficient! Maybe by stopping to use that vector of bool arrays...!
 */
void PermutableSpiketrainBuffer::clear()
{
	for (NeuronID ni = 0 ; ni < N_total ; ++ni)
	{
		for (AurynTime ti = 0 ; ti < timesteps_total ; ++ti)
		{
			theSpiketrains[ni][ti] = false;
		}
	}
	write_position = 0;
	read_position = 0;

	// reset permutation indices:
	if (permutationIndices.size() != N_total)
	{
		cout << "permutationIndices has changed!!" << endl;
	}
	permutationIndices.resize(N_total);  // should not be needed
	for (size_t i = 0; i != N_total; ++i)
		permutationIndices[i] = i;
}






void PermutableSpiketrainBuffer::permuteSpiketrains(SpikeContainer* whoToPermute, SpikeContainer* thePermutationIndices)
{
	if (whoToPermute->size() != thePermutationIndices->size())
	{
		cout << "Error: the sizes of whoToPermute and thePermutationIndices do not match!!" << endl;
	}


	for (unsigned int nii = 0 ; nii < whoToPermute->size() ; ++nii)
	{
		// assuming that thePermutationIndices is an exact permutation of whoToPermute:
		//cout << "Changing index " << (*whoToPermute)[nii] << " to " << (*thePermutationIndices)[nii] << endl;
		NeuronID theOldIndex = (*whoToPermute)[ nii ];
		NeuronID theNewIndex = (*whoToPermute)[ (*thePermutationIndices)[nii] ];
		//cout << "Changing index " << theOldIndex << " to " << theNewIndex << endl;
		//cout << "  size of permutationIndices is " << permutationIndices.size() << endl;
		//cout << "  the previous index at vector-field " << theOldIndex << " was: " << permutationIndices[theOldIndex] << endl;
		permutationIndices[theOldIndex] = theNewIndex;
		//permutationIndices[theNewIndex] = theOldIndex;
	}
	ispermuted = true;
}

void PermutableSpiketrainBuffer::getUnitsThatSpikeAtLeastOnce(SpikeContainer* targetspikecontainer, NeuronID N_presenting)
{
	for (NeuronID ni = 0 ; ni < N_presenting ; ++ni)
	{
		bool spikedAtLeastOnce = false;
		for (AurynTime ti = 0 ; ti < timesteps_total ; ++ti)
		{
			spikedAtLeastOnce = spikedAtLeastOnce || theSpiketrains[ni][ti];
		}
		if (spikedAtLeastOnce)
		{
			targetspikecontainer->push_back(ni);
		}
	}
}

//void PermutableSpiketrainBuffer::getActiveUnitsAndRandomSpiketimes(SpikeContainer* targetspikecontainer, vector<AurynTime>* randomspiketimes, NeuronID N_presenting, boost::random::mt19937 poly_gen, boost::random::uniform_int_distribution<> int_dist)
//{
//	for (NeuronID ni = 0 ; ni < N_presenting ; ++ni)
//	{
//		temp_spikepositions.clear();
//		bool spikedAtLeastOnce = false;
//		for (AurynTime ti = 0 ; ti < timesteps_total ; ++ti)
//		{
//			if (hasSpikeAt(ni,ti))
//			{
//				temp_spikepositions.push_back(ti);
//			}
//			spikedAtLeastOnce = spikedAtLeastOnce || theSpiketrains[ni][ti];  // faster than checking for vector.size later?
//		}
//		if (spikedAtLeastOnce)
//		{
//			targetspikecontainer->push_back(ni);
//
//			unsigned int latencyID = temp_spikepositions.size();
//			boost::random::uniform_int_distribution<>::param_type newrange = boost::random::uniform_int_distribution<>::param_type(0,latencyID);
//			int_dist.param(newrange);
//			AurynTime latency = temp_spikepositions[int_dist(poly_gen)];
//			randomspiketimes->push_back(latency);
//		}
//	}
//}

/**
 * Go through the pre-generated spiketrain of each named NeuronID, and choose a single random spike for each given neuron, and return its latency.
 * @param randomspiketimes The latencies of the single randomly chosen spike per given neuron.
 * @param activeUnits The neurons to check. The given neurons must be ensured to fire at least once. This can be ensured by using getUnitsThatSpikeAtLeastOnce() to generate inputs.
 * @param poly_gen The random number generator.
 * @param int_dist The int distribution to use. The range of the distribution will be adjusted for each spiketrain, because the number of spikes per spiketrain varies.
 */
void PermutableSpiketrainBuffer::getRandomSpiketimeOfActiveUnits(vector<AurynTime>* randomspiketimes, const SpikeContainer* activeUnits)
{
	//cout << " Entering getRandomSpiketimeOfActiveUnits().." << endl;
	for (auto nid : *activeUnits) // C++11 syntax.
	{
		temp_spikepositions.clear();
		for (AurynTime ti = 0 ; ti < timesteps_total ; ++ti)
		{
			if (theSpiketrains[nid][ti])
			{
				temp_spikepositions.push_back(ti); // store spike latency relative to beginning of buffer.
			}
		}
		// Now choosing a random spike within all spikes of this neuron, and writing the latency into the return vector:
		//   Performance todo: What's more costly: The random number generator or an "if" statement to exclude single-spike randomisation?
		unsigned int totalSpikesOfThisNeuronWithinBuffer = temp_spikepositions.size();
		unsigned int random_spike_id = 0;
		if (totalSpikesOfThisNeuronWithinBuffer > 1)
		{
			boost::random::uniform_int_distribution<int>::param_type newrange = boost::random::uniform_int_distribution<int>::param_type(0,totalSpikesOfThisNeuronWithinBuffer-1);
			int_dist.param(newrange);  // set range parameters
			random_spike_id = int_dist(temp_gen);
		}
		//else: random_spike_id remains at 0;


		//std::random_device rd;     // only used once to initialise (seed) engine
		//std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
		//std::uniform_int_distribution<int> uni(0,totalSpikesOfThisNeuronWithinBuffer-1); // guaranteed unbiased
		//auto random_spike_id = uni(rng);

		AurynTime latency = temp_spikepositions[random_spike_id];
		if (latency > timesteps_total)
		{
			// something went badly wrong!
			cout << "  Found a weird spike latency of " << latency << " timesteps! The random_spike_id was " << random_spike_id << " . The size of temp_spikepositions was " << temp_spikepositions.size() << " ."<< endl;
		}
		randomspiketimes->push_back(latency);
	}
}



















