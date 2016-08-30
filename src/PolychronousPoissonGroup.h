/*
 * PolychronousPoissonGroup.h
 *
 *  Created on: 07.04.2016
 *      Author: simon
 */

#ifndef POLYCHRONOUSPOISSONGROUP_H_
#define POLYCHRONOUSPOISSONGROUP_H_

#include "auryn/auryn_definitions.h"
#include "auryn/System.h"
#include "auryn/PoissonGroup.h"
#include "PermutableSpiketrainBuffer.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

using namespace std;

namespace auryn
{
//typedef vector< vector<bool> > PermutableSpiketrainBuffer; //!< Multiple SpikeContainers, one for each time step.
	typedef unsigned int PatternID; //!< NeuronID is an unsigned integeger type used to index neurons in Auryn.

/**
 * \brief A group that reshuffles Poisson noise spike trains to generate polychronous patterns.
 *
 * More specific description will follow.
 */
	class PolychronousPoissonGroup: public PoissonGroup
	{
	private:

		boost::random::mt19937 poly_gen;
		boost::random::uniform_real_distribution<> float_dist;

		vector<PermutableSpiketrainBuffer> buffers; // we want multiple buffers, so that one can be read from, one can be ready, and possibly a third could be receiving "live" data from another thread (in the future).
		unsigned int current_read_buffer;
		unsigned int current_write_buffer;
		vector<AurynFloat> representedValues;

		// here?
		bool issorted;
		bool isready;
		void sortByRank();
		void imprintPattern();

		string patterntimesfilename;
		ofstream patterntimesfile;

		AurynTime next_event;
		bool stimulus_active;
		int current_stimulus;
		vector<PatternID> stimuli_immediate;

		SpikeContainer spikingUnitIDs;
		vector<AurynTime> unorderedLatencies;

		void init(NeuronID N_presenting, NeuronID N_subpresenting,
				  AurynFloat duration, AurynFloat interval, NeuronID stimuli, string outputfilename);
		void initBuffers(int delaysteps);
		void generateNoiseAhead(PermutableSpiketrainBuffer* buffer, AurynTime steps);

		void checkAndUpdateTestingProtocol();

		vector<AurynFloat> testprotocolDurations;
		vector<AurynFloat> testprotocolPatternintervals;
		AurynTime next_phase_clock;
		int next_phase_id;


	protected:
		SpikeContainer pregen_spikes;  // this is used by distribute_spike() to trick PoissonGroup::evolve() into not sending the stuff upwards yet.
		//SpikeContainer evolve_spikes;  // this is used by PoissonGroup::evolve(). Separate from pregen_spikes for possible future multi-threading.
		virtual void distribute_spike(NeuronID theSpikingNeuron); ///< Allows the distribution target to be overwritten by sub-classes.
		virtual void constructNextBufferContents(PermutableSpiketrainBuffer* theBuffer);
		virtual void argsort(vector<AurynTime>* unorderedLatencies, SpikeContainer* orderingIndices);

	public:

		//NeuronID& N_total; // rename the size var for internal use
		NeuronID N_presenting;
		NeuronID N_subpresenting;

		bool twoLegged;	// should synchain-like patterns be one long chain (twoLegged=false) or more like a double-sided wave?
		bool useRandomPermutations; // or should patterns be random permutations of firing order?
		AurynFloat participationProbability; // to further hide patterns in background data
		NeuronID numPatterns;	// total number of distinct patterns (=stimuli)
		AurynTime patternDuration;	// may vary during runtime, if implemented
		AurynTime patternInterval;	// may vary during runtime, if implemented
		AurynTime max_patternDuration;	// fixed buffer size

		AurynTime max_patternInterval;	// fixed buffer size



		PolychronousPoissonGroup(NeuronID N_total, NeuronID N_presenting, NeuronID N_subpresenting,
								 AurynFloat duration, AurynFloat interval, NeuronID num_stimuli = 1,  AurynDouble rate=5. ,
								 string outputfilename = "stimulus.dat" );
		virtual ~PolychronousPoissonGroup();
		virtual void evolve();

		vector<PatternID> get_stimuli_immediate();

		virtual void seed(int s);

		PatternID getNumPatterns();
		AurynTime getMaxPatternDuration();
		void setMaxPatternDuration(AurynTime max_patternDuration);
		AurynTime getMaxPatternInterval() ;
		void setMaxPatternInterval(AurynTime max_patternInterval);

		template<typename T> class CompareIndicesByAnotherVectorValues
		{
			std::vector<T>* _values;
		public:
			CompareIndicesByAnotherVectorValues(std::vector<T>* values) : _values(values) {}
		public:
			bool operator()(const int& a, const int& b) const { return (*_values)[a] < (*_values)[b]; }
		};

		void setTestingProtocol(vector<AurynFloat> testprotocolDurations, vector<AurynFloat> testprotocolPatternintervals);

	};
}




#endif /* POLYCHRONOUSPOISSONGROUP_H_ */







