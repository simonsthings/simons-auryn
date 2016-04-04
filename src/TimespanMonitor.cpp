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

#include "TimespanMonitor.h"

TimespanMonitor::TimespanMonitor(string filename) : Monitor(filename)
{
}


TimespanMonitor::~TimespanMonitor()
{
	delete [] startclocks;
	delete [] stopclocks;
}



void TimespanMonitor::record_conditional()
{
	if (use_given_recordingtimes)
	{
		// record data only if within a predefined recording window (and the end of recording windows has not yet been reached)
		if ( (recordingtime_id < num_recording_timespans) && (sys->get_clock() >= startclocks[recordingtime_id]) )
		{
			// record data of this timestep:
			record_data();

			// advance to next recording window if indicated:
			if (sys->get_clock() >= stopclocks[recordingtime_id])
			{
				++recordingtime_id;
			}
		}
		// else do nothing (=do not record)
	}
	else
	{
		// just record everything.
		record_data();
	}
}


void TimespanMonitor::set_recording_times(auryn_vector_float* starttimes, auryn_vector_float* stoptimes)
{
	num_recording_timespans = starttimes->size;
	delete [] startclocks;  // don't care about NULL pointers. delete already takes care of that.
	delete [] stopclocks;
	startclocks = new unsigned int[num_recording_timespans];
	stopclocks = new unsigned int[num_recording_timespans];

	for (unsigned int i = 0; i < num_recording_timespans; ++i)
	{
		startclocks[i] = int(starttimes->data[i] / dt);  // convert from SI time in s to integer simulation steps
		stopclocks[i]  = int(stoptimes->data[i]  / dt);  // convert from SI time in s to integer simulation steps
	}

	use_given_recordingtimes = true;
	recordingtime_id = 0;
}


void TimespanMonitor::propagate()
{
	record_conditional();
}








