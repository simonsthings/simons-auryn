
#include <SpikeResponseMonitor.h>

SpikeResponseMonitor::SpikeResponseMonitor(PolychronousPoissonGroup *spikepatternprovider, NeuronGroup *responder,
										   SpikeContainer trackedResponseNeurons, string filename,
										   AurynTime binsize)
		: Monitor(filename)
{
	init(spikepatternprovider, responder, trackedResponseNeurons, filename, binsize);
}

SpikeResponseMonitor::~SpikeResponseMonitor()
{
}

void SpikeResponseMonitor::init(PolychronousPoissonGroup *theSpikepatternprovider, NeuronGroup *theResponder,
								SpikeContainer theTrackedResponseNeurons, string filename, AurynTime binsize)
{
	sys->register_monitor(this);

	spikepatternprovider = theSpikepatternprovider;
	responder = theResponder;
	maxPatternInterval = theSpikepatternprovider->getMaxPatternInterval();
	ssize = binsize;
	trackedResponseNeurons = theTrackedResponseNeurons;  // what does this actually do? Copy? Or just grasp on to this new object?

	PatternID numPatterns = theSpikepatternprovider->getNumPatterns();

	//responseTrackers.resize(theTrackedResponseNeurons.size());
	//for (int j = 0; j < maxNumPatterns; ++j)
	//for(auto iter : responseTrackers)
	//vector<std::vector<LatencyContainer> >::iterator iter;
	for (auto iter = theTrackedResponseNeurons.begin() ; iter != theTrackedResponseNeurons.end() ; ++iter)
	{
		// ensure each pattern's response tracker vectors are large enough to track each response neuron:
		//(*iter).resize(numPatterns);

		NeuronID theTrackedNeuronID = *iter;
		responseTrackers[theTrackedNeuronID]; // create new field. Throw away any return value here.
	}
	lastPatternResets.resize(numPatterns);

	outfile << setiosflags(ios::scientific) << setprecision(6);
}


void SpikeResponseMonitor::propagate()
{
	const vector<PatternID> &anyPatterns = spikepatternprovider->get_stimuli_immediate();
	//SpikeContainer* anyResponses = responder->get_spikes_immediate();
	const SpikeContainer &anyResponses = (*responder->get_spikes_immediate());

	// handle the arrival of a new pattern (by resetting the according timebin pointer).
	for (auto patternIter = anyPatterns.begin() ; patternIter != anyPatterns.end() ; ++patternIter)
	{
		//const PatternID thePatternID = *patternIter;
		lastPatternResets[*patternIter] = sys->get_clock();
	}

	auto totalTimesteps = 100.0 /dt;
	if (sys->get_clock() > totalTimesteps/2)
	{
		// put any new spikes into the appropriate bin:
		for (auto neuronIter = (anyResponses).begin() ; neuronIter != (anyResponses).end() ; ++neuronIter )
		{ // loop over incoming response neuronIDs:
			NeuronID theNeuronID = *neuronIter;
			auto& responseTrackersForThisNeuron = responseTrackers[theNeuronID];

			for (int pid = 0 ; pid < lastPatternResets.size() ; ++pid)
			{ // loop over all expected patterns:
				LatencyContainer &lc = responseTrackersForThisNeuron[pid];
				AurynTime latencyRelativeToThisPatternsOnset = sys->get_clock() - lastPatternResets[pid];
				if (latencyRelativeToThisPatternsOnset <= maxPatternInterval)
					lc[latencyRelativeToThisPatternsOnset]++;
				else
					// should never be reached if PolychronousPoissonGroup and this class have the same number of patterns set:
					cout << "A new pattern for pid " << pid << " has not appeared for longer than the given maxPatternInterval! Latency by now: " << latencyRelativeToThisPatternsOnset << " timesteps. maxPatternInterval: " << maxPatternInterval << ". Ignoring..." << endl;
			}
		}
	}


	if (sys->get_clock() % 10000 == 0)
	{
		displaySpikeCounts();
	}


}


void SpikeResponseMonitor::displaySpikeCounts()
{
	for (auto iter = responseTrackers.begin() ; iter != responseTrackers.end() ; ++iter)
	{
		auto& responseTrackersForThisNeuron = responseTrackers[*iter];
	}
}

//void SpikeResponseMonitor::resetTracker()
//{
//}
