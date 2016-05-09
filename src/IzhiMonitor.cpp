/* 
* Copyright 2014-2015 Friedemann Zenke
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
*
* If you are using Auryn or parts of it for your work please cite:
* Zenke, F. and Gerstner, W., 2014. Limits to high-speed simulations 
* of spiking neural networks using general-purpose computers. 
* Front Neuroinform 8, 76. doi: 10.3389/fninf.2014.00076
*/

#include "IzhiMonitor.h"

IzhiMonitor::IzhiMonitor(Izhikevich2003Group * source, NeuronID id, string filename, AurynTime stepsize) : TimespanMonitor(filename)
{
	init(source,id,filename,stepsize);
}

IzhiMonitor::~IzhiMonitor()
{
}

void IzhiMonitor::init(Izhikevich2003Group * source, NeuronID id, string filename, AurynTime stepsize)
{
	sys->register_monitor(this);

	src = source;
	ssize = stepsize;
	if ( ssize < 1 ) ssize = 1;

	nid = id;
	outfile << setiosflags(ios::scientific) << setprecision(6);
	outfile << "# time    ";
	outfile << "\t\t ampa   ";
	outfile << "\t\t nmda   ";
	outfile << "\t\t t_exc   ";
	outfile << "\t\t t_inh   ";
	outfile << "\t\t mem     ";
	outfile << "\t\t v_temp  ";
	outfile << "\t\t  u      ";
	outfile << "\t\t u_temp  ";
	outfile << "\t\ttempMem[0]";
	outfile << "\t\ttempMem[1]";
	outfile << "\t\ttempMem[2]";
	outfile << "\t\ttempMem[3]";
	outfile << "\t\ttempMem[4]";
	outfile << "\t\ttempMem[5]";
	outfile << "\t\ttempMem[6]";
	outfile << "\t\ttempMem[7]";
	outfile << "\t\ttempMem[8]";
	outfile << "\t\ttempMem[9]";
	outfile << "\t\ttempMem[10]";
	outfile << "\t\ttempMem[11]";
	outfile << "\t\ttempMem[12]";
	outfile << endl;
}

void IzhiMonitor::record_data()
{
	if (sys->get_clock()%ssize==0)
	{
		outfile << dt*(sys->get_clock()) << " \t" ;
		outfile << src->get_ampa(nid) << " \t" ;
		outfile << src->get_nmda(nid) << " \t" ;
		outfile << src->get_t_exc(nid) << " \t";
		outfile << src->get_t_inh(nid) << " \t";
		outfile << src->get_mem(nid) << " \t";
		outfile << src->get_v_temp(nid) << " \t";
		outfile << src->get_u(nid) << "\t";
		outfile << src->get_u_temp(nid) << "\t";
		outfile << src->get_tempMemState(0) << "\t";
		outfile << src->get_tempMemState(1) << "\t";
		outfile << src->get_tempMemState(2) << "\t";
		outfile << src->get_tempMemState(3) << "\t";
		outfile << src->get_tempMemState(4) << "\t";
		outfile << src->get_tempMemState(5) << "\t";
		outfile << src->get_tempMemState(6) << "\t";
		outfile << src->get_tempMemState(7) << "\t";
		outfile << src->get_tempMemState(8) << "\t";
		outfile << src->get_tempMemState(9) << "\t";
		outfile << src->get_tempMemState(10) << "\t";
		outfile << src->get_tempMemState(11) << "\t";
		outfile << src->get_tempMemState(12) << "\t";
		outfile << endl;
	}
}









