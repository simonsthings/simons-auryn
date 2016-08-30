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

#ifndef TIMESPANMONITOR_H_
#define TIMESPANMONITOR_H_

#include "auryn/auryn_definitions.h"
#include "auryn/Monitor.h"
#include "auryn/System.h"
#include "auryn/Connection.h"
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

namespace auryn
{
/*! \brief Adds a way to record only predefined timespans of data.
 *
 * If set_recording_times() has not been set (or you have overwritten propagate() and record_conditional() is not being used),
 * TimespanMonitor behaves just as Monitor. In your propagate function, call record_conditional() instead of directly writing
 * the output data. Implement your output data writing in the record_data() function.
 */

	class TimespanMonitor : protected Monitor
	{
	protected:
		unsigned int num_recording_timespans=0;
		unsigned int* startclocks = NULL;  // considering std::vector<uint> here. Votes?
		unsigned int* stopclocks  = NULL;  // considering std::vector<uint> here. Votes?
		bool use_given_recordingtimes = false;
		int recordingtime_id = 0;
		void record_conditional();
		virtual void record_data() = 0;
		/*! Standard constructor with file name*/
		TimespanMonitor(string filename);
		/*! Standard destructor  */
		virtual ~TimespanMonitor();

	public:
		/** Specify the times during which recording should be active. TODO: Just use normal std::vector<float> here. */
		void set_recording_times(auryn_vector_float* starttimes, auryn_vector_float* stoptimes);
		/*! Default implementation for propagate() */
		virtual void propagate();
	};
}



#endif /*TIMESPANMONITOR_H_*/
