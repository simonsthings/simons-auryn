
#include <StimulusResponseMonitor.h>

StimulusResponseMonitor::StimulusResponseMonitor(PolychronousPoissonGroup *spikepatternprovider, NeuronGroup *responder,
										   int trackFirstNeurons, string filename,
										   AurynTime binsize)
		: Monitor(filename)
{
	init(spikepatternprovider, responder, trackFirstNeurons, filename, binsize);
}

StimulusResponseMonitor::~StimulusResponseMonitor()
{
}

void StimulusResponseMonitor::init(PolychronousPoissonGroup *theSpikepatternprovider, NeuronGroup *theResponder,
								int trackFirstNeurons, string filename, AurynTime binsize)
{
	sys->register_monitor(this);

	spikepatternprovider = theSpikepatternprovider;
	responder = theResponder;
	maxPatternInterval = theSpikepatternprovider->getMaxPatternInterval();
	ssize = binsize;
	numTrackedResponseNeurons = trackFirstNeurons;  // what does this actually do? Copy? Or just grasp on to this new object?

	PatternID numPatterns = theSpikepatternprovider->getNumPatterns();
	requestedPatternPresentationsPerTrackingWindow = 50; // 50 = 10 seconds if there are 5 patterns per second.

	responseTrackers.resize(trackFirstNeurons);
	//for (int j = 0; j < maxNumPatterns; ++j)
	//for(auto iter : responseTrackers)
	vector<std::vector<LatencyContainer> >::iterator iter;
	for (iter = responseTrackers.begin() ; iter != responseTrackers.end() ; ++iter)
	{
		// ensure each pattern's response tracker vectors are large enough to track each response neuron:
		(*iter).resize(numPatterns);
	}
	lastPatternResets.resize(numPatterns);


	// statistics storage vars: //

	patternPresentationsInCurrentTrackingWindow.resize(numPatterns);

	singlePatternResponses.resize(trackFirstNeurons);
	for (auto & subvector : singlePatternResponses)
	{
		subvector.resize(numPatterns);
	}


	latencyDependentSpikeCountsInCurrentWindow.resize(trackFirstNeurons);
	for (auto & subvector : latencyDependentSpikeCountsInCurrentWindow)
	{
		subvector.resize(numPatterns);
	}

	truepositivesInCurrentTrackingWindow.resize(trackFirstNeurons);
	for (auto & subvector : truepositivesInCurrentTrackingWindow)
	{
		subvector.resize(numPatterns);
	}

	falsepositivesInCurrentTrackingWindow.resize(trackFirstNeurons);
	for (auto & subvector : falsepositivesInCurrentTrackingWindow)
	{
		subvector.resize(numPatterns);
	}


	outfile << "#  true_positive_rate  false_positive_rate  tpr/fpr " << endl;
	outfile << setiosflags(ios::scientific) << setprecision(6);
}


void StimulusResponseMonitor::propagate()
{
	const vector<PatternID> &anyPatterns = spikepatternprovider->get_stimuli_immediate();
	//SpikeContainer* anyResponses = responder->get_spikes_immediate();
	const SpikeContainer &anyResponses = (*responder->get_spikes_immediate());

	//if (anyPatterns.size() > 0 ) cout << "A Pattern was presented!" << endl;
	//if (anyResponses.size() > 0 ) cout << "At least one spike was fired!" << endl;

	// handle the arrival of a new pattern (by resetting the according timebin pointer).
	for (PatternID thePatternID : anyPatterns)
	{
		patternPresentationsInCurrentTrackingWindow[thePatternID]++;

		// TODO do the computations here! (well in some function of course) Also check if we should stop the evaluations now or continue for more evaluation windows.
		computePerPatternStatistics(thePatternID);

		if (patternPresentationsInCurrentTrackingWindow[thePatternID] > requestedPatternPresentationsPerTrackingWindow)
		{
			computeMultiPatternStatistics(thePatternID);
			resetMultiPatternData(thePatternID);
		}



		// reset last pattern occurrence time:
		lastPatternResets[thePatternID] = sys->get_clock();

		// reset responses for new pattern that is just starting:
		for ( auto neuronIter = singlePatternResponses.begin() ; neuronIter != singlePatternResponses.end() ; neuronIter++  )
		{
			auto spikeTimesForEachPattern = (*neuronIter);
			vector<AurynTime> theSpikeTimes = spikeTimesForEachPattern[thePatternID];
			theSpikeTimes.clear();
		}

	}

	auto totalTimesteps = 100.0 /dt;
	if (sys->get_clock() > totalTimesteps/2)
	{
		// put any new spikes into the appropriate bin:
		for (auto neuronIter = (anyResponses).begin() ; neuronIter != (anyResponses).end() ; ++neuronIter )
		{ // loop over incoming response neuronIDs:
			NeuronID theNeuronID = *neuronIter;
			auto& responseTrackersForThisNeuron = responseTrackers[*neuronIter];

			for (int pid = 0 ; pid < lastPatternResets.size() ; ++pid)
			{ // loop over all expected patterns:
				LatencyContainer &lc = responseTrackersForThisNeuron[pid];
				AurynTime latencyRelativeToThisPatternsOnset = sys->get_clock() - lastPatternResets[pid];
				if (latencyRelativeToThisPatternsOnset <= maxPatternInterval)
				{
					// record multi-pattern peak data:
					lc[latencyRelativeToThisPatternsOnset]++;

					// record single-pattern spike response data:
					singlePatternResponses[theNeuronID][pid].push_back(latencyRelativeToThisPatternsOnset);

					// if latency is within X ms after pattern onset, do:
					//truePositives[theNeuronID][pid]++;

					// if latency is longer than X ms after pattern onset, do:
					//falsePositives[theNeuronID][pid]++; // DANGER: multiple spikes possible!! If not compensated, this would lead to wrong results!
				}
				else
					// should never be reached if PolychronousPoissonGroup and this class have the same number of patterns set:
					cout << "A new pattern for pid " << pid << " has not appeared for longer than the given maxPatternInterval! Latency by now: " << latencyRelativeToThisPatternsOnset << " timesteps. maxPatternInterval: " << maxPatternInterval << ". Ignoring..." << endl;
			}
		}
	}


	if (sys->get_clock() % long(25.0/dt) == 0)
	{
		displaySpikeCounts();
	}


}


void StimulusResponseMonitor::displaySpikeCounts()
{
	//cout << "Displaying spike counts!" << endl;
	for (auto iter = responseTrackers.begin() ; iter != responseTrackers.end() ; ++iter)
	{
		auto & responseTrackersForThisNeuron = *iter;
		for (int pid = 0 ; pid < lastPatternResets.size() ; ++pid)
		{
			LatencyContainer &lc = responseTrackersForThisNeuron[pid];
			unordered_map<unsigned int, unsigned int>::iterator latencyIter;
			//cout << "Spike response delays (nX, pX): [";

			for (int ti = 0; ti < maxPatternInterval; ++ti)
			{
				//cout << lc[ti] << " ";
			}
			//cout << "]" << endl;

			lc.clear();
		}
	}
}

void StimulusResponseMonitor::computePerPatternStatistics(PatternID thePatternID)
{
	// TODO: use data on single-pattern responses to find out typical number of responses and false-positive rates etc.

	AurynTime givenTPmaxtime = (AurynTime)  (50.0 *1e-3/dt); // e.g. 50 ms/dt
	AurynTime givenFPmintime = (AurynTime) (150.0 *1e-3/dt); // e.g. 150 ms/dt

//	vector<AurynTime> latencyDependentSpikeCountLatencies;
//	latencyDependentSpikeCountLatencies.push_back( (AurynTime) 10.0 *1e-3/dt );
//	latencyDependentSpikeCountLatencies.push_back( (AurynTime) 20.0 *1e-3/dt );
//	latencyDependentSpikeCountLatencies.push_back( (AurynTime) 30.0 *1e-3/dt );
//	latencyDependentSpikeCountLatencies.push_back( (AurynTime) 40.0 *1e-3/dt );


	// for each response neuron:
	for (int ni = 0; ni < singlePatternResponses.size(); ++ni)
	{
		auto spikeTimesForEachPattern = singlePatternResponses[ni];
		vector<AurynTime> theSpikeTimes = spikeTimesForEachPattern[thePatternID];


//		vector<SpikeCount> & latencyDependentSpikeCount = latencyDependentSpikeCountsInCurrentWindow[ni][thePatternID];
//		latencyDependentSpikeCount.clear();
//		latencyDependentSpikeCount.resize(latencyDependentSpikeCountLatencies.size());

		//spikeCountInCurrentWindow[ni][thePatternID].push_back(theSpikeTimes.size());

		//cout << "Spiketimes (length: " <<  theSpikeTimes.size()  << "): ";
		bool truepositive = false;
		bool falsepositive = false;
		for (AurynTime theSpiketime : theSpikeTimes)
		{
			// track true positives and false positives:
			if (theSpiketime < givenTPmaxtime) truepositive = true;
			if (theSpiketime > givenFPmintime) falsepositive = true;

			//cout << theSpiketime << " " << endl;
//			// track latency-dependent firing rates:
//			bool previouslyAssigned = false; // works like an else if, but compatible with the loop below!
//			for (int li = 0; li < latencyDependentSpikeCountLatencies.size(); ++li)
//				if ( (theSpiketime < latencyDependentSpikeCountLatencies[li]) && (!previouslyAssigned) )
//					latencyDependentSpikeCount[li]++;
		}
		//cout << endl;
		if (truepositive)   truepositivesInCurrentTrackingWindow[ni][thePatternID]++;
		if (falsepositive) falsepositivesInCurrentTrackingWindow[ni][thePatternID]++;

		// reset single pattern responses because we have now dealt with them.
		singlePatternResponses[ni][thePatternID].clear();
	}

}


void StimulusResponseMonitor::computeMultiPatternStatistics(PatternID thePatternID)
{
	// TODO: use the responseTrackers to find peaks of spike responses within the last X pattern presentations!

	const unsigned int &numPatternPresentations = patternPresentationsInCurrentTrackingWindow[thePatternID];

	for (int ni = 0; ni < truepositivesInCurrentTrackingWindow.size(); ++ni)
	{
		//cout << "numPatternPresentations: " << numPatternPresentations;
		//cout << ", TP: " <<  truepositivesInCurrentTrackingWindow[ni][thePatternID];
		//bv cout << ", FP: " << falsepositivesInCurrentTrackingWindow[ni][thePatternID];

		float tpr = ( (float)truepositivesInCurrentTrackingWindow[ni][thePatternID] / (float)numPatternPresentations );
		float fpr = ( (float)falsepositivesInCurrentTrackingWindow[ni][thePatternID] / (float)numPatternPresentations );

		//cout << setprecision(4);
		//cout << "  True positive rate: " << tpr << ", " << "False positive rate: " << fpr << ". " << endl;

		outfile << tpr << " "  << fpr << " "  << tpr/fpr << " " << endl;
	}

}


void StimulusResponseMonitor::resetMultiPatternData(PatternID thePatternID)
{
	patternPresentationsInCurrentTrackingWindow[thePatternID] = 0;

	// clear temp data (hope this does not hurt performence too much)
	for (int ni = 0; ni < truepositivesInCurrentTrackingWindow.size(); ++ni)
	{
		truepositivesInCurrentTrackingWindow[ni][thePatternID] = 0;
		falsepositivesInCurrentTrackingWindow[ni][thePatternID] = 0;
	}
}


















