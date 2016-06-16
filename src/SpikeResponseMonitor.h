
#ifndef SPIKERESPONSEMONITOR_H_
#define SPIKERESPONSEMONITOR_H_

#include "auryn_definitions.h"
#include "Monitor.h"
#include "TimespanMonitor.h"
#include "PolychronousPoissonGroup.h"
#include "System.h"
#include "Connection.h"
#include <fstream>
#include <iomanip>
#include <unordered_map>

using namespace std;

typedef unsigned int SpikeCount;
typedef unordered_map<AurynTime,SpikeCount> LatencyContainer;

/*! \brief Records the response delays and other metrics of one specific unit
 * from a NeuronGroup in response to pattern presentation by a PolychronousPoissonGroup.
 *
 * Metrics recorded shall include:
 * - (mean) number of response spikes per pattern presentation over the last X timebins
 * - median delay for 1st, 2nd, 3rd, etc response spike
 * - speed of delay change for each peak (as a measure of response timing stability)
 * - mean firing rates of the responding neuron during eary and late pattern presentation, and during noise intervals
 * - response reliability:
 * -- yes/no to pattern presentation
 * --
 **/
class SpikeResponseMonitor : public Monitor
{
protected:
	int numTrackedResponseNeurons; // delete this?
	PolychronousPoissonGroup* spikepatternprovider;
	NeuronGroup* responder;
	AurynTime maxPatternInterval;
	AurynTime ssize;

	vector<vector<LatencyContainer> > responseTrackers;
	vector<vector<unsigned int>::iterator> responseTrackerIterators;
	vector<AurynTime> lastPatternResets;

	vector<vector<vector<AurynTime> > > singlePatternResponses;  // a SpikeContainer for each stimulus/pattern type.

	vector<vector<vector<SpikeCount> > > latencyDependentSpikeCountsInCurrentWindow;
	vector<PatternID> patternPresentationsInCurrentTrackingWindow;
	vector<vector<unsigned int> > truepositivesInCurrentTrackingWindow;
	vector<vector<unsigned int> > falsepositivesInCurrentTrackingWindow;

	void init(PolychronousPoissonGroup *spikepatternprovider, NeuronGroup *responder,
			  int trackedResponseNeurons, string filename, AurynTime binsize);

public:
	SpikeResponseMonitor(PolychronousPoissonGroup *spikepatternprovider, NeuronGroup *responder,
						 int trackFirstNeurons, string filename, AurynTime binsize);
	virtual ~SpikeResponseMonitor();

	virtual void propagate() override;
	void displaySpikeCounts();

	void computePeakStatistics(PatternID thePatternID);

	void computeSpikeStatistics(PatternID thePatternID);

};

#endif /*SPIKERESPONSEMONITOR_H_*/
