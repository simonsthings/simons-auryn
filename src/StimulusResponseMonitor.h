
#ifndef SPIKERESPONSEMONITOR_H_
#define SPIKERESPONSEMONITOR_H_

#include "auryn/auryn_definitions.h"
#include "auryn/Monitor.h"
#include "TimespanMonitor.h"
#include "PolychronousPoissonGroup.h"
#include "auryn/System.h"
#include "auryn/Connection.h"
#include <fstream>
#include <iomanip>
#include <unordered_map>

using namespace std;

namespace auryn
{
	typedef unsigned int SpikeCount;
	typedef unsigned int PatternCount;
	typedef unordered_map<AurynTime,SpikeCount> LatencyContainer;

/*! \brief Records the response delays and other metrics of one specific unit
 * from a NeuronGroup in response to pattern presentation by a PolychronousPoissonGroup.
 *
 * Metrics recorded shall include:
 * - (mean) number of response spikes per pattern presentation over the last X timebins
 * - median delay for 1st, 2nd, 3rd, etc response spike
 * - speed of delay change for each peak (as a measure of response timing stability)
 * - mean firing rates of the responding neuron during early and late pattern presentation, and during noise intervals
 * - response reliability:
 * -- yes/no to pattern presentation
 * --
 **/
	class StimulusResponseMonitor : public Monitor
	{
	protected:
		PolychronousPoissonGroup* spikepatternprovider;		//!< Inputs
		NeuronGroup* responder;								//!< Outputs

		NeuronID neuronToTrack;		//!< Which Neuron to track with this instance of StimulusResponseMonitor.
		PatternID patternToTrack;	//!< Which Pattern to track with this instance of StimulusResponseMonitor.

		LatencyContainer responseTracker;  //!< Used to track all lags of all responses relative to each pattern onset. Will be used to compute response peaks.
		AurynTime lastPatternArrivalClock;

		vector<AurynTime> singlePatternResponses;  //!< This tracks the latencies of all responses between two pattern presentations.

		PatternCount requestedPatternPresentationsPerTrackingWindow;
		PatternCount patternPresentationsInCurrentTrackingWindow;
		unsigned int truepositivesInCurrentTrackingWindow;
		unsigned int falsepositivesInCurrentTrackingWindow;

		void init(PolychronousPoissonGroup *spikepatternprovider, NeuronGroup *responder,
				  NeuronID neuronToTrack, PatternID patternToTrack, PatternID binsizeInPatternPresentations);

	public:
		StimulusResponseMonitor(PolychronousPoissonGroup *spikepatternprovider, NeuronGroup *responder, string filename, AurynTime binsizeInPatternPresentations,
								NeuronID neuronToTrack, PatternID patternToTrack);
		virtual ~StimulusResponseMonitor();

		virtual void propagate() override;

		void writeTFPNRatesToFile();

		/** For each pattern presentation, this computes true positives and false positives. The result is then stored in the fields truepositivesInCurrentTrackingWindow and falsepositivesInCurrentTrackingWindow. */
		void computeHitsAndFalsealarms();

		void resetTFPNs();
	};

}


#endif /*SPIKERESPONSEMONITOR_H_*/
