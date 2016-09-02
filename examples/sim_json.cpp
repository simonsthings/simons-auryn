/* 
*
* This file is based on the Auryn example files, a simulation package for plastic
* spiking neural networks.
* 
* Auryn is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* Auryn is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with Auryn.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "auryn.h"
#include "../src/auryn/System.h"
#include <ctime>
#include <boost/property_tree/json_parser.hpp>
#include <StimulusResponseMonitor.h>

using namespace std;
using namespace auryn;

namespace po = boost::program_options;
namespace mpi = boost::mpi;

void defineDefaultParameters(boost::property_tree::ptree const& pt);
void print(boost::property_tree::ptree const& pt, string indent );
float getSamplinginterval(boost::property_tree::ptree const& pt, string path);
void readParametersFromJSONfileAndCMDline(int ac, char* av[], boost::property_tree::ptree & simparams);
//void setupHistoryTracking(SpikingGroup* poisson, NeuronGroup* detector_neuron, STDPConnection* con1, const boost::property_tree::ptree& simparams);
SpikingGroup* setupPresynapticGroup(const boost::property_tree::ptree& simparams);
NeuronGroup* setupPostsynapticGroup(const boost::property_tree::ptree& simparams);
DuplexConnection* setupConnection(SpikingGroup* poisson, NeuronGroup* detector_neuron, const boost::property_tree::ptree& simparams);
void defineGlobals(mpi::communicator world, mpi::environment& env, const boost::property_tree::ptree& simparams);
int main(int ac, char* av[]);



void defineDefaultParameters(boost::property_tree::ptree const& pt)
{
}

void print(boost::property_tree::ptree const& pt, string indent = "  ")
{
    using boost::property_tree::ptree;
    ptree::const_iterator end = pt.end();
    for (ptree::const_iterator it = pt.begin(); it != end; ++it) {
        std::cout << indent << it->first << ": " << it->second.get_value<std::string>() << std::endl;
        print(it->second,indent+"  ");
    }
}

float getSamplinginterval(boost::property_tree::ptree const& pt, string path)
{
	if (pt.get<string>(path) == "dt")
	{
		return (float) auryn_timestep;
	}
	else
	{
		return pt.get<float>(path);
	}
}

void readParametersFromJSONfileAndCMDline(int ac, char* av[], boost::property_tree::ptree & simparams)
{
	string settingsfile_simulation = "settings_simulation.json";

	// BEGIN Command line parameters
    try {

        po::options_description desc("Allowed options");
        desc.add_options()
			("help", "produce help message")
            ("settingsfile", po::value<string>(), "alternate simulation settings filename and path. (default: ./settings_simulation.json)")
			("show-settings", "display simulation settings as read from json file")
			("write-example-settings", "write example settings to two files called 'settings_simulation_example.json'")
            ("general.simtime", po::value<double>(), "override the simulation time")
            ("general.outfileprefix", po::value<string>(), "override the output file prefix")
            //("Npre", po::value<int>(), "number of input units (= poisson group size)")
            //("seed", po::value<int>(), "(example) overwrite random seed")
        ;


        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << desc << "\n";
            exit(1);
        }

        if (vm.count("write-example-settings")) {
            cout << "Writing example settings to './settings_simulation_example.json' and exiting.";
            cout << " (TODO: actually do this!)";
            cout << "\n";
            exit(1);
        }

        if (vm.count("settingsfile")) {
            std::cout << "Reading simulation settings from file: "
                 << vm["settingsfile"].as<string>() << ".\n";
            settingsfile_simulation = vm["settingsfile"].as<string>();
        }

        // BEGIN testing json parser
        boost::property_tree::ptree pt;  // the property tree for storing the simulation parameters!
    	std::ifstream jsonfile(settingsfile_simulation);
    	if (jsonfile)
    	{
    		std::stringstream ss;
    		ss << jsonfile.rdbuf();
    		try
    		{
    			boost::property_tree::read_json(ss, pt);
    			cout << "Successfully parsed json file." << endl;
    			simparams.swap(pt);
    		}
    		catch (std::exception const& e)
    		{
    			std::cerr << "The error is: " << e.what() << std::endl;
    			cout << "Using default parameters because something went wrong with the json file parsing." << endl;
    		}
    	}
    	else
    	{
    		cout << "Warning: could not open json file \"" << settingsfile_simulation << "\" for some reason. (TODO: make this a proper log message!)" << endl;
    	}
        // END testing json parser

        if (vm.count("show-settings")) {
            cout << "Using these simulation settings:" << "\n";
//            showsettings = true;
    		//cout << "Now showing data from freshly parsed json file:" << endl;
    		print(simparams);
    		//cout << "Finished showing that data!" << endl;
        }

//        if (vm.count("prefix")) {
//            std::cout << "file prefix set to "
//                 << vm["prefix"].as<string>() << ".\n";
//            file_prefix = vm["prefix"].as<string>();
//        }
//
//        if (vm.count("simtime")) {
//            cout << "simtime set to "
//                 << vm["simtime"].as<double>() << " seconds.\n";
//			simtime = vm["simtime"].as<double>();
//        }

    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        std::exit(1);
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        std::exit(1);
    }
	// END Command line parameters

}

template <typename T>
std::vector<T> as_vector(boost::property_tree::ptree const& pt, boost::property_tree::ptree::key_type const& key)
{
	std::vector<T> r;
	for (auto& item : pt.get_child(key))
		r.push_back(item.second.get_value<T>());
	return r;
}

SpikingGroup* setupPresynapticGroup(const boost::property_tree::ptree& simparams)
{
	SpikingGroup* poisson;

	try
	{

		NeuronID Npre = simparams.get<NeuronID>("neurongroups.inputs.N");
		string requestedNeuronGroupClass = simparams.get<string>("neurongroups.inputs.type");
		cout << "The requested input group is: " << requestedNeuronGroupClass << endl;
		if ("PoissonGroup" == requestedNeuronGroupClass)
		{
			AurynDouble inputpoprate = simparams.get<double>("neurongroups.inputs.rate");
			poisson = new PoissonGroup(Npre, inputpoprate);
			((PoissonGroup*) (poisson))->seed(simparams.get<unsigned int>("neurongroups.inputs.randomseed"));
		}
		else if ("FileInputGroup" == requestedNeuronGroupClass)
		{
			string inputRASfile = simparams.get<string>("neurongroups.inputs.rasfilename");
			cout << "The input RAS file: " << inputRASfile << endl;
			poisson = new FileInputGroup(Npre, inputRASfile,false,0.0);
		}
		else if ("StructuredPoissonGroup" == requestedNeuronGroupClass)
		{
			AurynFloat patDuration = simparams.get<AurynFloat>("neurongroups.inputs.patternduration");
			AurynFloat patInterval = simparams.get<AurynFloat>("neurongroups.inputs.patterninterval");
			unsigned int numStimuli = simparams.get<unsigned int>("neurongroups.inputs.numberofstimuli");
			AurynDouble inputpoprate = simparams.get<double>("neurongroups.inputs.rate");
			string patOccurrencesFilename = simparams.get<string>("neurongroups.inputs.patternOccurrencesFilename");
			poisson = new StructuredPoissonGroup(Npre, patDuration, patInterval, numStimuli, inputpoprate, patOccurrencesFilename);
			((StructuredPoissonGroup*) (poisson))->seed(simparams.get<unsigned int>("neurongroups.inputs.randomseed"));
		}
		else if ("PolychronousPoissonGroup" == requestedNeuronGroupClass)
		{
			NeuronID Npre_presenting = simparams.get<NeuronID>("neurongroups.inputs.N_presenting");
			NeuronID Npre_subpresenting = simparams.get<NeuronID>("neurongroups.inputs.N_subpresenting");
			AurynFloat patDuration = simparams.get<AurynFloat>("neurongroups.inputs.patternduration");
			AurynFloat patInterval = simparams.get<AurynFloat>("neurongroups.inputs.patterninterval");
			unsigned int numStimuli = simparams.get<unsigned int>("neurongroups.inputs.numberofstimuli");
			AurynDouble inputpoprate = simparams.get<double>("neurongroups.inputs.rate");
			string patOccurrencesFilename = simparams.get<string>("neurongroups.inputs.patternOccurrencesFilename");
			poisson = new PolychronousPoissonGroup(Npre, Npre_presenting, Npre_subpresenting, patDuration, patInterval, numStimuli, inputpoprate,
					patOccurrencesFilename);
			((PolychronousPoissonGroup*) (poisson))->seed(simparams.get<int>("neurongroups.inputs.randomseed"));


			for (auto i : as_vector<AurynFloat>(simparams, "general.testingProtocol.durations")) std::cout << i << ' ';
			std::cout << '\n';
			for (auto j : as_vector<AurynFloat>(simparams, "general.testingProtocol.intervals")) std::cout << j << ' ';


			vector<AurynFloat> testprotocolDurations = as_vector<AurynFloat>(simparams, "general.testingProtocol.durations");
			//testprotocolDurations.push_back(30);
			//testprotocolDurations.push_back(30);
			//testprotocolDurations.push_back(simparams.get<AurynFloat>("general.simtime")-60);
			vector<AurynFloat> testprotocolPatternintervals = as_vector<AurynFloat>(simparams, "general.testingProtocol.intervals");
			//testprotocolPatternintervals.push_back(0.3);   // negative values indicate noise only.
			//testprotocolPatternintervals.push_back(0.2);  // in seconds.
			//testprotocolPatternintervals.push_back(10.0); // in seconds.
			((PolychronousPoissonGroup*) (poisson))->setTestingProtocol(testprotocolDurations,testprotocolPatternintervals);
		}
		else
		{
			throw AurynGenericException();
		}

	}

    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        std::exit(1);
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        std::exit(1);
    }


	return poisson;
}

NeuronGroup* setupPostsynapticGroup(const boost::property_tree::ptree& simparams)
{
	NeuronGroup* postGroup;
	try
	{


		NeuronID N_post = simparams.get<NeuronID>("neurongroups.outputs.N");
		Izhikevich2003Group* detector_neuron = new Izhikevich2003Group(N_post);
		detector_neuron->set_projMult( simparams.get<float>("neurongroups.outputs.projMult") );
		detector_neuron->use_recovery = simparams.get<bool>("neurongroups.outputs.userecovery");

		cout << "IzhikevichGroup.projMult: " << detector_neuron->get_projMult() << endl;
		cout << "IzhikevichGroup.use_recovery: " << detector_neuron->use_recovery << endl;
		postGroup = detector_neuron;


		//	IFGroup* detector_neuron = new IFGroup(N_post);
		//	cout << "IFNeuron.ampa: " << detector_neuron->get_ampa(0) << endl;
		//	cout << "IFNeuron.cursyn: " << detector_neuron->get_cursyn(0) << endl;
		//	cout << "IFNeuron.effective_load: " << detector_neuron->get_effective_load() << endl;
		//	cout << "IFNeuron.gaba: " << detector_neuron->get_gaba(0) << endl;
		//	cout << "IFNeuron.locked_range: " << detector_neuron->get_locked_range() << endl;
		//	cout << "IFNeuron.locked_rank: " << detector_neuron->get_locked_rank() << endl;
		//	cout << "IFNeuron.mem: " << detector_neuron->get_mem(0) << endl;pych
		//	cout << "IFNeuron.name: " << detector_neuron->get_name() << endl;
		//	cout << "IFNeuron.nmda: " << detector_neuron->get_nmda(0) << endl;
		//	cout << "IFNeuron.num_spike_attributes: " << detector_neuron->get_num_spike_attributes() << endl;
		//	cout << "IFNeuron.post_size: " << detector_neuron->get_post_size() << endl;
		//	cout << "IFNeuron.pre_size: " << detector_neuron->get_pre_size() << endl;
		//	cout << "IFNeuron.rank_size: " << detector_neuron->get_rank_size() << endl;
		//	cout << "IFNeuron.size: " << detector_neuron->get_size() << endl;
		//	cout << "IFNeuron.tau_ampa: " << detector_neuron->get_tau_ampa() << endl;
		//	cout << "IFNeuron.tau_gaba: " << detector_neuron->get_tau_gaba() << endl;
		//	cout << "IFNeuron.tau_mem: " << detector_neuron->get_tau_mem() << endl;
		//	cout << "IFNeuron.tau_nmda: " << detector_neuron->get_tau_nmda() << endl;
		//	cout << "IFNeuron.uid: " << detector_neuron->get_uid() << endl;
		//	cout << "IFNeuron.vec_size: " << detector_neuron->get_vector_size() << endl;
		//
		//	cout << "Adjusting default values of IFNeuron..." << endl;
		//	detector_neuron->set_tau_mem(0.4);
		//
		//	cout << "IFNeuron.tau_mem: " << detector_neuron->get_tau_mem() << endl;


		//	CubaIFGroup* detector_neuron = new CubaIFGroup(N_post);  // doesn't even start. Need more starting weight??
		//	cout << "CubaIFGroup.ampa: " << detector_neuron->get_ampa(0) << endl;
		//	cout << "CubaIFGroup.cursyn: " << detector_neuron->get_cursyn(0) << endl;
		//	cout << "CubaIFGroup.effective_load: " << detector_neuron->get_effective_load() << endl;
		//	cout << "CubaIFGroup.gaba: " << detector_neuron->get_gaba(0) << endl;
		//	cout << "CubaIFGroup.locked_range: " << detector_neuron->get_locked_range() << endl;
		//	cout << "CubaIFGroup.locked_rank: " << detector_neuron->get_locked_rank() << endl;
		//	cout << "CubaIFGroup.mem: " << detector_neuron->get_mem(0) << endl;
		//	cout << "CubaIFGroup.name: " << detector_neuron->get_name() << endl;
		//	cout << "CubaIFGroup.nmda: " << detector_neuron->get_nmda(0) << endl;
		//	cout << "CubaIFGroup.num_spike_attributes: " << detector_neuron->get_num_spike_attributes() << endl;
		//	cout << "CubaIFGroup.post_size: " << detector_neuron->get_post_size() << endl;
		//	cout << "CubaIFGroup.pre_size: " << detector_neuron->get_pre_size() << endl;
		//	cout << "CubaIFGroup.rank_size: " << detector_neuron->get_rank_size() << endl;
		//	cout << "CubaIFGroup.size: " << detector_neuron->get_size() << endl;
		//	cout << "CubaIFGroup.uid: " << detector_neuron->get_uid() << endl;
		//	cout << "CubaIFGroup.vec_size: " << detector_neuron->get_vector_size() << endl;
		//
		//	cout << "CubaIFGroup.bg_current: " << detector_neuron->get_bg_current(0) << endl;
		//
		//	cout << "Adjusting default values of IFNeuron..." << endl;
		//	detector_neuron->set_all_bg_currents(0.011);
		//
		//	cout << "CubaIFGroup.bg_current: " << detector_neuron->get_bg_current(0) << endl;

			//SIFGroup* detector_neuron = new SIFGroup(N_post);   // type not found? It is included everywhere?
			//TIFGroup* detector_neuron = new TIFGroup(N_post);
			//IafPscDeltaGroup* detector_neuron = new IafPscDeltaGroup(N_post);
			//AdExGroup* detector_neuron = new AdExGroup(N_post);
			//AIFGroup* detector_neuron = new AIFGroup(N_post);
			//AIF2Group* detector_neuron = new AIF2Group(N_post);   // this is the only Group that stops! But I guess some other process is involved here...


	}
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        std::exit(1);
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        std::exit(1);
    }

    return postGroup;

}

DuplexConnection* setupConnection(SpikingGroup* poisson, NeuronGroup* detector_neuron, const boost::property_tree::ptree& simparams)
{
	DuplexConnection* con1;

	try
	{

		// END neuron groups & settings.
		// BEGIN Connections & settings:
		AurynWeight weight = simparams.get<float>("connectionsets.con1.initialweight"); //0.1;
		//AurynFloat sparseness = 0.99999999;
		AurynFloat sparseness = 1.0;
		AurynFloat learningrate = simparams.get<float>("connectionsets.con1.stdprule.learningrate"); //0.01;  // learningrate (?)
		AurynFloat A_plus = simparams.get<float>("connectionsets.con1.stdprule.A_plus");
		AurynFloat A_minus = simparams.get<float>("connectionsets.con1.stdprule.A_minus");
		AurynFloat tau_pre = simparams.get<float>("connectionsets.con1.stdprule.tau_plus");
		AurynFloat tau_post = simparams.get<float>("connectionsets.con1.stdprule.tau_minus");
		AurynFloat maxweight = simparams.get<float>("connectionsets.con1.maximumweight"); //1.0f;
		TransmitterType transmitter = GLUT;

		int STDPruleID = 0;

		switch (STDPruleID)
		{
			case 0:
			{
				STDPConnection* theConn = new STDPConnection(poisson, detector_neuron, weight, sparseness, learningrate, tau_pre, tau_post, maxweight, transmitter, "testSTDP");
				double ms = 1e-3;  // milisecond scale. Need to find a pretty way of handling this.
				theConn->A *= A_minus;// *dt/ms;
				theConn->B *= A_plus;// *dt/ms;
				con1 = theConn;
			};break;
			case 1:
			{
				STDPwdConnection* theConn = new STDPwdConnection(poisson, detector_neuron, weight, sparseness, learningrate, maxweight, transmitter, "testSTDP");
				//theConn->set_alphalambda_alternative(A_plus,A_minus,learningrate);
				theConn->set_mu_plus(1.0);
				theConn->set_mu_minus(1.0);
				con1 = theConn;
			};break;
			case 2:
			{
//				STDPlsConnection* theConn = new STDPlsConnection(poisson, detector_neuron, weight, sparseness, learningrate, maxweight, transmitter, "testSTDP");
				//theConn->set_alphalambda_alternative(A_plus,A_minus,learningrate);
				//theConn->set_mu_plus(0.0);
				//theConn->set_mu_minus(0.0);
//				con1 = theConn;
			};break;
			default:
			{
				throw std::invalid_argument( "received invalid STDPruleID." );
			};break;
		}

	}
	catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		std::exit(1);
	}
	catch(...) {
		cerr << "Exception of unknown type!\n";
		std::exit(1);
	}

	return con1;
}

void setupDetailedHistoryTracking(SpikingGroup *poisson, NeuronGroup *detector_neuron,
								  DuplexConnection *con1, const boost::property_tree::ptree &simparams)
{
	string tmpstr;


	// END Connections & settings.
	// BEGIN History setup (writing to disk):
	auryn_vector_float* starttimes = auryn_vector_float_alloc(100);
	auryn_vector_float* stoptimes = auryn_vector_float_alloc(100);
	stringstream ssStart(simparams.get<string>("recordings.dtintervalsAsStrings.starttimes"));
	stringstream ssStop(simparams.get<string>("recordings.dtintervalsAsStrings.stoptimes"));
	string buf1, buf2;
	unsigned int counter = 0;
	while (ssStart >> buf1)
	{
		ssStop >> buf2;
		starttimes->data[counter] = (AurynFloat) (stof(buf1));
		stoptimes->data[counter] = (AurynFloat) (stof(buf2));
		++counter;
		cout << "Reading in starttime: " << buf1 << " ...and Stoptime: " << buf2
		<< endl;
	}
	starttimes->size = counter;
	stoptimes->size = counter;

	// TODO: actually use this instead of the parsed strings above:
	vector<AurynFloat> theStarttimes = as_vector<AurynFloat>(simparams, "recordings.dtintervalsAsFloats.starttimes");
	vector<AurynFloat> theStoptimes  = as_vector<AurynFloat>(simparams, "recordings.dtintervalsAsFloats.stoptimes");


	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += "_inputs.ras";
	SpikeMonitor* smon_in = new SpikeMonitor(poisson, tmpstr.c_str());
	smon_in->set_recording_times(starttimes, stoptimes);

	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += "_inputs.spk";
	BinarySpikeMonitor* binsmon_in = new BinarySpikeMonitor(poisson, tmpstr.c_str());

	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += "_inputs.poprate";
	PopulationRateMonitor* rmon_in = new PopulationRateMonitor(poisson,
															   tmpstr.c_str(),
															   getSamplinginterval(simparams,
																				   "recordings.inputs.samplinginterval_poprate"));





	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += "_outputs.ras";
	SpikeMonitor* smon_out = new SpikeMonitor(detector_neuron, tmpstr.c_str());
	smon_out->set_recording_times(starttimes, stoptimes);

	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += "_outputs.spk";
	BinarySpikeMonitor* binsmon_out = new BinarySpikeMonitor(detector_neuron, tmpstr.c_str());

	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += "_outputs.rate";
	RateMonitor* rmon_out = new RateMonitor(detector_neuron, tmpstr.c_str(),
											getSamplinginterval(simparams,
																"recordings.outputs.samplinginterval_rate"));


	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += ".mem";
	VoltageMonitor* vmon = new VoltageMonitor(detector_neuron, 0, tmpstr.c_str(), getSamplinginterval(simparams, "recordings.outputs.samplinginterval_membranes"));
	vmon->set_recording_times(starttimes, stoptimes);

	/*
	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += ".ampa";
	AmpaMonitor* amon = new AmpaMonitor(detector_neuron, 0, tmpstr.c_str(),
				(AurynTime)(getSamplinginterval(simparams, "recordings.outputs.samplinginterval_ampa") / (float)dt)); // divide by dt because constructor interfaces are not (yet) normed across classes.
	amon->set_recording_times(starttimes, stoptimes);


	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += ".nmda";
	NmdaMonitor* nmon = new NmdaMonitor(detector_neuron, 0, tmpstr.c_str(),
				(AurynTime)(getSamplinginterval(simparams, "recordings.outputs.samplinginterval_nmda") / dt)); // divide by dt because constructor interfaces are not (yet) normed across classes.
	nmon->set_recording_times(starttimes, stoptimes);
	*/




	//tmpstr = simparams.get<string>("general.outfileprefix");
	//tmpstr += ".testing";
	//TestingMonitor * tmon = new TestingMonitor( tmpstr.c_str() ); // divide by dt because constructor interfaces are not (yet) normed across classes.
	//tmon->set_recording_times(starttimes,stoptimes);


	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += ".izhi";
	IzhiMonitor* izhmon = new IzhiMonitor((Izhikevich2003Group*)detector_neuron, 0, tmpstr.c_str(),
				  (AurynTime)(getSamplinginterval(simparams,"recordings.outputs.samplinginterval_ampa") / auryn_timestep)); // divide by dt because constructor interfaces are not (yet) normed across classes.
	izhmon->set_recording_times(starttimes, stoptimes);






	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += ".weightsum";
	WeightSumMonitor* weightsum = new WeightSumMonitor(con1, tmpstr.c_str(),
													   getSamplinginterval(simparams,
																		   "recordings.con1.samplinginterval_weightsum"));


	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += ".weightstats";
	WeightStatsMonitor* weightstats = new WeightStatsMonitor(con1,
															 tmpstr.c_str(),
															 getSamplinginterval(simparams,
																				 "recordings.con1.samplinginterval_weightstats"));


	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += ".weightmatrix";
	WeightMatrixMonitor* weightmatrix = new WeightMatrixMonitor(con1,tmpstr.c_str(),
																getSamplinginterval(simparams,"recordings.con1.samplinginterval_weightmatrix"));


	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += ".stimulusdetectionstatistics.txt";
	int numTrackedNeurons = 1;
	int binsize_patternPresentations = 50;
	StimulusResponseMonitor* srm = new StimulusResponseMonitor((PolychronousPoissonGroup *) poisson, detector_neuron, tmpstr, binsize_patternPresentations, 0, 0);



}

void setupReducedHistoryTracking(SpikingGroup *poisson, NeuronGroup *detector_neuron, DuplexConnection *con1, const boost::property_tree::ptree &simparams)
{
	string tmpstr;

	tmpstr = simparams.get<string>("general.outfileprefix");
	tmpstr += ".stimulusdetectionstatistics.txt";
	int numTrackedNeurons = 1;
	StimulusResponseMonitor* srm = new StimulusResponseMonitor((PolychronousPoissonGroup *) poisson, detector_neuron, tmpstr, 50, 0, 0);
}


void setupHistoryTracking(SpikingGroup *poisson, NeuronGroup *detector_neuron, DuplexConnection *con1, const boost::property_tree::ptree &simparams)
{
	bool useDetailedTracking; // system will exit on error, so no default declaration needed.

	try {useDetailedTracking = simparams.get<bool>("recordings.detailedtracking");}
	catch(exception& e) { cerr << "error: " << e.what() << "\n"; std::exit(1); }
	catch(...) { cerr << "Exception of unknown type!\n"; std::exit(1); }

	if (useDetailedTracking)
		setupDetailedHistoryTracking(poisson, detector_neuron, con1, simparams);
	else
		setupReducedHistoryTracking(poisson, detector_neuron, con1, simparams);
}

void defineGlobals(mpi::communicator world, mpi::environment& env, const boost::property_tree::ptree& simparams)
{
	mpicommunicator = &world;
	try
	{
		std::ostringstream out;
		out << simparams.get<string>("general.outfileprefix") << "." << world.rank() << ".log";
		logger = new Logger(out.str(), world.rank(), PROGRESS, EVERYTHING);
	}
	catch (AurynOpenFileException& excpt)
	{
		std::cerr << "Cannot proceed without log file. Exiting all ranks ..." << '\n';
		env.abort(1);
	}
	catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		std::exit(1);
	}
	catch(...) {
		cerr << "Exception of unknown type!\n";
		std::exit(1);
	}
	sys = new System(&world);
}

int main(int ac, char* av[])
{
	int errcode = 0;

	//bool showsettings = false;
	//char strbuf [255];
	//std::istringstream settingsfile_recordings;
    boost::property_tree::ptree simparams;  // the property tree for storing the simulation parameters!
    /* Fill params structure with default parameters! These will also be printed to an example file if requested. */
    defineDefaultParameters(simparams);
    readParametersFromJSONfileAndCMDline(ac,av,simparams);


	// BEGIN Global definitions
	//mpi::environment env(ac, av);
	//mpi::communicator world;
	//defineGlobals(world, env, simparams);
	auryn_init(ac, av);
	// END Global definitions



	// BEGIN neuron groups & settings:
	SpikingGroup* poisson = setupPresynapticGroup(simparams);
	NeuronGroup* detector_neuron = setupPostsynapticGroup(simparams);
	// END neuron groups & settings.



	// BEGIN Connections & settings:
	DuplexConnection* con1 = setupConnection(poisson, detector_neuron, simparams);
	// END Connections & settings.



	// BEGIN History setup (writing to disk):
	setupHistoryTracking(poisson, detector_neuron, con1, simparams);
	// END History setup (writing to disk).



	// BEGIN running the simulation:
	logger->msg("Running ...",PROGRESS);
	sys->run(simparams.get<float>("general.simtime"));
	// END running the simulation.


	if (errcode)
		mpienv->abort(errcode);
	logger->msg("Freeing ...",PROGRESS,true);
	auryn_free();
	return errcode;
}

