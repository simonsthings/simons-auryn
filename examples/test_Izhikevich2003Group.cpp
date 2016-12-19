/* 
* Copyright 2014 Friedemann Zenke
*
* This file is part of Auryn, a simulation package for plastic
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

using namespace std;
using namespace auryn;

namespace po = boost::program_options;
namespace mpi = boost::mpi;

int main(int ac, char* av[]) 
{

	string dir = "./";
	string file_prefix = "test_Izhi2007";

	char strbuf [255];
	string msg;
	double mV = 1e-3;  // 1 / 1000 th.
	double pA = 1;  // need to check what Auryn actually expects here.

	NeuronID size = 1000;
	NeuronID seed = 1;
	double kappa = 5.;
	double simtime = 2.0;  // in seconds
	//double stepAt  = simtime/4.0;  // in seconds
	double stepAt  = 1.01;  // in seconds

	if (stepAt >= simtime)
		throw logic_error("the step time must be within the sim time!!");

	int errcode = 0;

    try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("simtime", po::value<double>(), "simulation time")
            ("kappa", po::value<double>(), "poisson group rate")
            ("dir", po::value<string>(), "output directory")
            ("size", po::value<int>(), "poisson group size")
            ("seed", po::value<int>(), "random seed")
        ;

        po::variables_map vm;        
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }

        if (vm.count("kappa")) {
            cout << "kappa set to " 
                 << vm["kappa"].as<double>() << ".\n";
			kappa = vm["kappa"].as<double>();
        } 

        if (vm.count("dir")) {
            std::cout << "dir set to " 
                 << vm["dir"].as<string>() << ".\n";
			dir = vm["dir"].as<string>();
        } 

        if (vm.count("simtime")) {
            cout << "simtime set to " 
                 << vm["simtime"].as<double>() << ".\n";
			simtime = vm["simtime"].as<double>();
        } 

        if (vm.count("size")) {
            cout << "size set to " 
                 << vm["size"].as<int>() << ".\n";
			size = vm["size"].as<int>();
        } 

        if (vm.count("seed")) {
            cout << "seed set to " 
                 << vm["seed"].as<int>() << ".\n";
			seed = vm["seed"].as<int>();
        } 
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

	// BEGIN Global stuff
	auryn_init(ac, av);
	// END Global stuff

	NeuronID N_post = 1;

	enum IzhiTypes  {Zenke_Izhi, Izhi2003, Izhi2007, IzhiAlternative};

	CurrentInjector* theInjector;

	//IzhiTypes theType = Izhi2003;
	IzhiTypes theType = IzhiAlternative;
	//IzhiTypes theType = Zenke_Izhi;

	if (theType == Izhi2003)
	{
		Izhikevich2003Group* izhi_neuron = new Izhikevich2003Group(N_post);
		//theInjector  = new CurrentInjector(izhi_neuron,"mem");
		//theInjector  = new CurrentInjector(izhi_neuron,"t_exc");
		theInjector  = new CurrentInjector(izhi_neuron,"g_ampa");
		theInjector->set_current(0,0.0);  // start with no added current


		VoltageMonitor* vmon = new VoltageMonitor(izhi_neuron, 0, "test_Izhi.mem.txt");
		StateMonitor* umon = new StateMonitor(izhi_neuron, 0,"u", "test_Izhi.u.txt");
		StateMonitor* exc_mon = new StateMonitor(izhi_neuron, 0,"t_exc", "test_Izhi.t_exc.txt");
		SpikeMonitor * smon_e = new SpikeMonitor( izhi_neuron, "test_Izhi.ras.txt", size);
	}
	else if (theType == IzhiAlternative)
	{
		IzhikevichAlternativeGroup* izhi_neuron = new IzhikevichAlternativeGroup(N_post);
		//theInjector  = new CurrentInjector(izhi_neuron,"mem");
		//theInjector  = new CurrentInjector(izhi_neuron,"t_exc");
		theInjector  = new CurrentInjector(izhi_neuron,"g_ampa");
		theInjector->set_current(0,0.0);  // start with no added current


		VoltageMonitor* vmon = new VoltageMonitor(izhi_neuron, 0, "test_Izhi.mem.txt");
		StateMonitor* umon = new StateMonitor(izhi_neuron, 0,"u", "test_Izhi.u.txt");
		StateMonitor* exc_mon = new StateMonitor(izhi_neuron, 0,"t_exc", "test_Izhi.t_exc.txt");
		SpikeMonitor * smon_e = new SpikeMonitor( izhi_neuron, "test_Izhi.ras.txt", size);
	}
	else if (theType == Izhi2007)
	{
		// do something here too, please.
	}
	else if (theType == Zenke_Izhi)
	{
		IzhikevichGroup* izhi_neuron = new IzhikevichGroup(N_post);

		izhi_neuron->avar = 0.02;    // adaptation variable rate constant
		izhi_neuron->bvar = 0.2;     // subtreshold adaptation
		izhi_neuron->cvar = -65e-3f; // reset voltage
		izhi_neuron->dvar = 6.0e-3;
		izhi_neuron->set_tau_ampa(1e-9);
		izhi_neuron->set_tau_gaba(1e-9);


		//theInjector  = new CurrentInjector(izhi_neuron,"mem");
		//theInjector  = new CurrentInjector(izhi_neuron,"t_exc");
		theInjector  = new CurrentInjector(izhi_neuron,"g_ampa");
		theInjector->set_current(0,0.0);  // start with no added current


		VoltageMonitor* vmon = new VoltageMonitor(izhi_neuron, 0, "test_Izhi.mem.txt");
		StateMonitor* umon = new StateMonitor(izhi_neuron, 0,"izhi_adaptation", "test_Izhi.u.txt");
		StateMonitor* exc_mon = new StateMonitor(izhi_neuron, 0,"cur_exc", "test_Izhi.t_exc.txt");
		StateMonitor* ampa_mon = new StateMonitor(izhi_neuron, 0,"g_ampa", "test_Izhi.g_ampa.txt");
		SpikeMonitor * smon_e = new SpikeMonitor( izhi_neuron, "test_Izhi.ras.txt", size);
	}

	///< simulate until the current step:
	if (!sys->run(stepAt,false))
			errcode = 2;


	double stepcurrent = 14.0 * 1e+3 ;//*pA; // pA=1 currently.  // still not sure which order of mag this needs to be...
	//double stepcurrent = 2.355 * 1e+3 ;//*pA; // pA=1 currently.  // still not sure which order of mag this needs to be...
	//double stepcurrent = 0.0;
	theInjector->set_current(0,stepcurrent);


	///< sim until the end:
	if (!sys->run(simtime-stepAt,false))
			errcode = 1;

	if (errcode)
		mpienv->abort(errcode);
	logger->msg("Freeing ...",PROGRESS,true);
	auryn_free();
	return errcode;
}













