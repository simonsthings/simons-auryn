
#include <StimulusResponseMonitor.h>

using namespace auryn;

StimulusResponseMonitor::StimulusResponseMonitor(PolychronousPoissonGroup *spikepatternprovider, NeuronGroup *responder, string filename,
												 AurynTime binsizeInPatternPresentations,
												 NeuronID neuronToTrack, PatternID patternToTrack)
		: Monitor(filename)
{
	init(spikepatternprovider, responder, neuronToTrack, patternToTrack, binsizeInPatternPresentations);
}

StimulusResponseMonitor::~StimulusResponseMonitor()
{
}

void StimulusResponseMonitor::init(PolychronousPoissonGroup *theSpikepatternprovider, NeuronGroup *theResponder,
								   NeuronID theNeuronToTrack, PatternID thePatternToTrack, PatternID binsizeInPatternPresentations)
{
	sys->register_device(this);

	spikepatternprovider = theSpikepatternprovider;
	responder = theResponder;
	//double binsizeInSeconds = binsizeInSeconds * dt;

	neuronToTrack = theNeuronToTrack;
	patternToTrack = thePatternToTrack;

	// Todo: use binsizeInClocksteps (=ssize) to compute this:
	//auto vfds = spikepatternprovider->patternInterval;
	requestedPatternPresentationsPerTrackingWindow = binsizeInPatternPresentations; // 50 = 10 seconds if there are 5 patterns per second.


	singlePatternResponses.clear();

	resetTFPNs();


	outfile << "# time  true_positive_rate  false_positive_rate  tpr-fpr  tpr/fpr " << endl;
	outfile << setiosflags(ios::fixed) << setprecision(6);
}


void StimulusResponseMonitor::propagate()
{
	// handle the arrival of a new pattern (by resetting the according timebin pointer).
	const vector<PatternID> &anyPatterns = spikepatternprovider->get_stimuli_immediate();
	for (PatternID thePatternID : anyPatterns)
	{
		if (thePatternID == patternToTrack)
		{
			patternPresentationsInCurrentTrackingWindow++;

			// TODO Check if we should stop the evaluations now or continue for more evaluation windows.
			computeHitsAndFalsealarms();

			if (patternPresentationsInCurrentTrackingWindow >= requestedPatternPresentationsPerTrackingWindow)
			{
				writeTFPNRatesToFile();
				resetTFPNs();
			}

			// reset last pattern occurrence time:
			lastPatternArrivalClock = sys->get_clock();

			// reset responses for new pattern that is just starting:
			singlePatternResponses.clear(); // I specifically only reset this counter if exactly this pattern is presented again.
		}
	}


	// handle the arrival of a new (response?) spike:
	const SpikeContainer &anyResponses = (*responder->get_spikes_immediate());
	for (NeuronID theNeuronID : anyResponses)
	{
		if ( (theNeuronID == neuronToTrack) && (lastPatternArrivalClock > 0) ) // check if correct neuron and previous patterns have arrived
		{
			// put any new spikes into the appropriate bin:
			AurynTime latencyRelativeToThisPatternsOnset = sys->get_clock() - lastPatternArrivalClock;
			if (latencyRelativeToThisPatternsOnset <= spikepatternprovider->getMaxPatternInterval())
			{
				// record multi-pattern peak data:
				responseTracker[latencyRelativeToThisPatternsOnset]++;

				// record single-pattern spike response data:
				singlePatternResponses.push_back(latencyRelativeToThisPatternsOnset);

				// if latency is within X ms after pattern onset, do:
				//truePositives[theNeuronID][pid]++;

				// if latency is longer than X ms after pattern onset, do:
				//falsePositives[theNeuronID][pid]++; // DANGER: multiple spikes possible!! If not compensated, this would lead to wrong results!
			}
			else
			{
				// should never be reached if PolychronousPoissonGroup and this class have the same number of patterns set:
				cout << "A new pattern for pid " << patternToTrack << " has not appeared for longer than the given maxPatternInterval! Latency by now: " << latencyRelativeToThisPatternsOnset << " timesteps. maxPatternInterval: " << spikepatternprovider->getMaxPatternInterval() << ". Ignoring..." << endl;
			}
		}
	}
}



void StimulusResponseMonitor::computeHitsAndFalsealarms()
{
	// TODO: use data on single-pattern responses to find out typical number of responses and false-positive rates etc.

	AurynTime givenTPmaxClock = (AurynTime)  (50.0 *1e-3/auryn_timestep); // e.g. 50 ms/dt
	AurynTime givenFPminClock = (AurynTime) (150.0 *1e-3/auryn_timestep); // e.g. 150 ms/dt

//	vector<AurynTime> latencyDependentSpikeCountLatencies;
//	latencyDependentSpikeCountLatencies.push_back( (AurynTime) 10.0 *1e-3/dt );
//	latencyDependentSpikeCountLatencies.push_back( (AurynTime) 20.0 *1e-3/dt );
//	latencyDependentSpikeCountLatencies.push_back( (AurynTime) 30.0 *1e-3/dt );
//	latencyDependentSpikeCountLatencies.push_back( (AurynTime) 40.0 *1e-3/dt );

//	vector<SpikeCount> & latencyDependentSpikeCount = latencyDependentSpikeCountsInCurrentWindow[ni][thePatternID];
//	latencyDependentSpikeCount.clear();
//	latencyDependentSpikeCount.resize(latencyDependentSpikeCountLatencies.size());

	//spikeCountInCurrentWindow[ni][thePatternID].push_back(theSpikeTimes.size());

	//cout << "Spiketimes (length: " <<  theSpikeTimes.size()  << "): ";
	bool truepositive = false;
	bool falsepositive = false;
	for (AurynTime theSpiketime : singlePatternResponses)
	{
		// track true positives and false positives:
		if (theSpiketime < givenTPmaxClock) truepositive = true;
		if (theSpiketime > givenFPminClock) falsepositive = true;

		//cout << theSpiketime << " " << endl;
//		// track latency-dependent firing rates:
//		bool previouslyAssigned = false; // works like an else if, but compatible with the loop below!
//		for (int li = 0; li < latencyDependentSpikeCountLatencies.size(); ++li)
//			if ( (theSpiketime < latencyDependentSpikeCountLatencies[li]) && (!previouslyAssigned) )
//				latencyDependentSpikeCount[li]++;
	}
	//cout << endl;
	if (truepositive)   truepositivesInCurrentTrackingWindow++;
	if (falsepositive) falsepositivesInCurrentTrackingWindow++;

	// reset single pattern responses because we have now dealt with them.
	singlePatternResponses.clear();


}

/**
 * When the preset number of patterns has been received, this function writes the resulting true positives rate to a file.
 */
void StimulusResponseMonitor::writeTFPNRatesToFile()
{
	// TODO: use the responseTrackers to find peaks of spike responses within the last X pattern presentations!

	//	const unsigned int &numPatternPresentations = patternPresentationsInCurrentTrackingWindow;

	//cout << "numPatternPresentations: " << numPatternPresentations;
	//cout << ", TP: " <<  truepositivesInCurrentTrackingWindow[ni][thePatternID];
	//bv cout << ", FP: " << falsepositivesInCurrentTrackingWindow[ni][thePatternID];

	float tpr = ( (float)truepositivesInCurrentTrackingWindow / (float)patternPresentationsInCurrentTrackingWindow );
	float fpr = ( (float)falsepositivesInCurrentTrackingWindow / (float)patternPresentationsInCurrentTrackingWindow );

	//cout << setprecision(4);
	//cout << "  True positive rate: " << tpr << ", " << "False positive rate: " << fpr << ". " << endl;

	outfile << sys->get_time() << " " << tpr << " "  << fpr << " "  << tpr-fpr << " "  << tpr/fpr << " " << endl;

}


void StimulusResponseMonitor::resetTFPNs()
{
	patternPresentationsInCurrentTrackingWindow = 0;

	truepositivesInCurrentTrackingWindow = 0;
	falsepositivesInCurrentTrackingWindow = 0;
}


















