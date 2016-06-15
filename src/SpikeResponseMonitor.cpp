
#include <SpikeResponseMonitor.h>

SpikeResponseMonitor::SpikeResponseMonitor(PolychronousPoissonGroup *spikepatternprovider, NeuronGroup *responder, SpikeContainer trackedResponseNeurons,
										   unsigned int maxPatternInterval, string filename, AurynTime binsize,
										   unsigned int maxNumPatterns)
		: Monitor(filename)
{
	init(spikepatternprovider, responder, trackedResponseNeurons, maxPatternInterval, filename, binsize, maxNumPatterns);
}

SpikeResponseMonitor::~SpikeResponseMonitor()
{
}

void SpikeResponseMonitor::init(PolychronousPoissonGroup *theSpikepatternprovider, NeuronGroup *theResponder,
								SpikeContainer theTrackedResponseNeurons, unsigned int theMaxPatternInterval, string filename, AurynTime binsize,
								unsigned int maxNumPatterns)
{
	sys->register_monitor(this);

	spikepatternprovider = theSpikepatternprovider;
	responder = theResponder;
	maxPatternInterval = theMaxPatternInterval;
	ssize = binsize;
	trackedResponseNeurons = theTrackedResponseNeurons;  // what does this actually do? Copy? Or just grasp on to this new object?

	responseTrackers.resize(theTrackedResponseNeurons.size());
	//for (int j = 0; j < maxNumPatterns; ++j)
	//for(auto iter : responseTrackers)
	vector<std::vector<LatencyContainer> >::iterator iter;
	for (iter = responseTrackers.begin() ; iter != responseTrackers.end() ; ++iter)
	{
		vector<LatencyContainer> &responseTrackersForThisNeuron = (*iter);
		//cout << "The size of the vector of latency containers is: " << responseTrackersForThisPattern.size() << " before resizing, ";
		// ensure each pattern's response tracker vectors are enough to track each response neuron:
		responseTrackersForThisNeuron.resize(maxNumPatterns);
		//cout << "and " << responseTrackersForThisPattern.size() << " afterwards." << endl;
	}
	lastPatternResets.resize(maxNumPatterns);

//
//	for (auto iterator = responseTrackers.begin() ; iterator != responseTrackers.end() ; iterator++ )
//	{
//		vector<LatencyContainer> &responseTrackersForOnePattern = *iterator;
//		for (int i = 0; i < trackedResponseNeurons.size(); ++i)
//		{
//			responseTrackersForOnePattern[i] = new
//		}
//		responseTrackerForOnePattern.resize(maxPatternInterval);
//		responseTrackerIterators.push_back(responseTrackerForOnePattern.begin());
//	}

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
//		vector<unsigned int> &responseTrackerForOnePattern = responseTrackers[*patternIter];
		//finaliseResponsesSinceLastPresentation(responseTrackerForOnePattern);
		//responseTrackerIterators[*patternIter] = responseTrackerForOnePattern.begin();
		//responseTrackerIterators[*patternIter]++; //

		//const PatternID thePatternID = *patternIter;
		lastPatternResets[*patternIter] = sys->get_clock();

//		for (auto neuronIter = (*anyResponses).begin() ; neuronIter != (*anyResponses).end() ; ++neuronIter )
//		{
//		}
	}


	// put any new spikes into the appropriate bin:
	for (auto neuronIter = (anyResponses).begin() ; neuronIter != (anyResponses).end() ; ++neuronIter )
	{ // loop over incoming response neuronIDs:
		NeuronID theNeuronID = *neuronIter;
		auto& responseTrackersForThisNeuron = responseTrackers[*neuronIter];

		for (int pid = 0; pid < lastPatternResets.size(); ++pid)
		{ // loop over all expected patterns:
			LatencyContainer &lc = responseTrackersForThisNeuron[pid];
			AurynTime latencyRelativeToThisPatternsOnset = sys->get_clock() - lastPatternResets[pid];
			if (latencyRelativeToThisPatternsOnset <= maxPatternInterval)
				lc[latencyRelativeToThisPatternsOnset]++;
			else
				cout << "A new pattern for pid " << pid << " has not appeared for longer than the given maxPatternInterval! Latency by now: " << latencyRelativeToThisPatternsOnset << " timesteps. maxPatternInterval: " << maxPatternInterval << ". Ignoring..." << endl;
		}
	}





	
	if (sys->get_clock()%ssize==0) {
		outfile << dt*(sys->get_clock()) << " " << responder->get_nmda(trackedResponseNeurons[0]) << "\n";
	}
}



