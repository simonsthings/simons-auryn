#!/usr/bin/python
import numpy as np
import pylab as pl
import struct


current_version = (0,8,0)

class AurynBinaryFile:
    '''
    This class is the abstract base class to access binary Auryn files.
    '''
    def unpack(self, data):
        return struct.unpack(self.data_format, data)

    def get_frame(self,idx):
        '''
        Returns decoded contents of a frame with index idx as tuple 
        '''
        pos = idx*self.frame_size
        self.datafile.seek(pos,0)
        data = self.datafile.read(self.frame_size)
        frame = self.unpack(data)
        return frame

    def find_frame(self, time=0, lower=False):
        '''
        Find frame index of given time by bisection.
        '''
        idx_lo = 1
        idx_hi = self.last_frame
        while idx_lo+1<idx_hi:
            pivot = (idx_lo+idx_hi)/2
            at, nid = self.get_frame(pivot)
            ft = at*self.timestep
            if ft>time:
                idx_hi = pivot
            else:
                idx_lo = pivot
        if lower:
            return idx_lo
        else:
            return idx_hi

    def read_frames(self, bufsize):
        return self.datafile.read(bufsize*self.frame_size)

    def refresh(self):
        ''' Reload file internal properties in case the file has changed'''

        # Reader spk file header
        data = self.datafile.read(self.frame_size)
        self.header = struct.unpack(self.data_format, data)
        self.timestep = 1.0/self.header[0]
        version_code = int(self.header[1])%1000
        # TODO add header signature checks
        self.file_version = (version_code%10, (version_code%100)/10, (version_code%1000)/100)
        if self.class_version != self.file_version:
            print("Warning! Version mismatch between the decoding tool and the file version.")
            print("AurynBinarySpikeFile %s"%str(self.class_version))
            print("Fileversion %s"%str(self.file_version))

        # Determine size of file
        self.datafile.seek(0,2)
        self.filesize = self.datafile.tell()
        self.last_frame = self.filesize/self.frame_size-1
        self.num_data_frames = self.last_frame-1
        self.datafile.seek(self.frame_size,0)
        if self.debug:
            print "Filesize", self.filesize
            print "Number of frames", self.num_data_frames

        # Determine min and max time
        at,val = self.get_frame(1)
        self.t_min = at*self.timestep
        at,val = self.get_frame(self.last_frame)
        self.t_max = at*self.timestep

    def open_file(self):

        # TODO add exception handling
        self.frame_size = struct.calcsize(self.data_format)
        self.datafile = open(self.filename, "rb")

        self.refresh()


class AurynBinaryStateFile(AurynBinaryFile):
    '''
    This class gives abstract access to binary Auryn state monitor file.

    Public methods:
    get_data: extracts time series data from the specified interval 
    '''
    def __init__(self, filename, debug_output=False):
        # These params might have to adapted ot different Auryn datatypes and versions
        self.data_format = "@If"
        self.class_version = current_version 

        self.debug = debug_output
        self.filename = filename
        self.open_file()

    def __del__(self):
        self.datafile.close()

    def get_data(self, t_start=0.0, t_stop=1e32):
        ''' Returns timeseries of state for given temporal interval'''
        idx_start = self.find_frame( t_start, lower=False )
        idx_stop  = self.find_frame( t_stop, lower=True )
        start_pos = idx_start*self.frame_size
        num_elements = idx_stop-idx_start

        self.datafile.seek(start_pos,0)
        raw_data = self.datafile.read(num_elements*self.frame_size)

        data = []
        for i in xrange(num_elements):
            at, val = struct.unpack_from(self.data_format, raw_data, i*self.frame_size)
            data.append((self.timestep*at, val))
        return data


class AurynBinarySpikeFile(AurynBinaryFile):
    '''
    This class gives abstract access to binary Auryn spike raster file (spk).

    Public methods:
    get_spikes: extracts spikes (tuples of time and neuron id) for a given temporal range.
    get_spike_times_from_interval: extracts the spike times of a single unit and a given temporal range.
    '''
    def __init__(self, filename, debug_output=False):
        # These params might have to adapted ot different Auryn datatypes and versions
        self.data_format = "@II"
        self.class_version = current_version

        self.debug = debug_output
        self.filename = filename
        self.open_file()

    def __del__(self):
        self.datafile.close()

    def get_spikes( self, t_start=0.0, t_stop=1e32, max_id=1e32 ):
        idx_start = self.find_frame( t_start, lower=False )
        idx_stop = self.find_frame( t_stop, lower=True )
        start_pos = idx_start*self.frame_size
        num_elements = idx_stop-idx_start

        self.datafile.seek(start_pos,0)
        data = self.datafile.read(num_elements*self.frame_size)

        spikes = []
        for i in xrange(num_elements):
            at, nid = struct.unpack_from(self.data_format, data, i*self.frame_size)
            if nid<max_id:
                spikes.append((self.timestep*at, nid))
        return spikes

    def get_spike_counts( self, t_start=0.0, t_stop=1e32, min_size=1 ):
        idx_start = self.find_frame( t_start, lower=False )
        idx_stop = self.find_frame( t_stop, lower=True )
        start_pos = idx_start*self.frame_size
        num_elements = idx_stop-idx_start

        self.datafile.seek(start_pos,0)
        data = self.datafile.read(num_elements*self.frame_size)

        counts = np.zeros(min_size)
        for i in xrange(num_elements):
            at, nid = struct.unpack_from(self.data_format, data, i*self.frame_size)
            if nid >= len(counts):
                counts.resize(nid+1)
            counts[nid] += 1
        return counts

    def get_firing_rates( self, t_start=0.0, t_stop=1e32, min_size=1 ):
        return self.get_spike_counts( t_start, t_stop, min_size )/(t_stop-t_start)

    def get_last( self, seconds=1.0 ):
        ''' Returns the last x seconds of spikes'''
        return self.get_spikes(t_start=self.t_max-seconds)

    def get_spike_times( self, neuron_id=0, t_start=0, t_stop=1e32 ):
        idx_start = self.find_frame( t_start, lower=False )
        idx_stop = self.find_frame( t_stop, lower=True )

        idx_cur = idx_start
        start_pos = idx_start*self.frame_size
        self.datafile.seek(start_pos,0)

        bufsize = 1024*1024
        spike_times = []
        while idx_cur < idx_stop:
            if idx_stop-idx_cur<bufsize:
                bufsize = idx_stop-idx_cur
            idx_cur = idx_cur+bufsize
            data = self.read_frames(bufsize)
            for i in xrange(bufsize):
                at, nid = struct.unpack_from(self.data_format, data, i*self.frame_size)
                if nid==neuron_id:
                    spike_times.append(self.timestep*at)
        return spike_times

class AurynBinarySpikeView:
    '''
    A wrapper class for easy extraction of spikes from multiple spk files from different ranks.
    '''
    def __init__(self, filenames):
        self.filenames = filenames
        self.spike_files = []
        for filename in self.filenames:
            self.spike_files.append( AurynBinarySpikeFile(filename) )

    def sort_spikes(self, spikes):
        spikes.sort(key=lambda tup: tup[0])

    def get_spikes( self, t_start=0.0, t_stop=1e32, max_id=1e32 ):
        spikes = []
        for spk in self.spike_files:
            spikes.extend( spk.get_spikes( t_start, t_stop, max_id ) )

        self.sort_spikes(spikes)
        return spikes

    def get_last( self, seconds=1.0 ):
        spikes = []
        for spk in self.spike_files:
            spikes.extend( spk.get_last( seconds ) )

        self.sort_spikes(spikes)
        return spikes

    def time_triggered_histogram(self, trigger_times, time_offset=0.0, time_window=100e-3, max_neuron_id=1024):
        '''
        Sums spikes within a given time window around given trigger times.

        This function can be used to compute STAs, reverese correlations, etc.. 
        for instance to compute a linear receptive field. 

        Keyword arguments:
        trigger_times -- list of trigger times
        time_window -- size of time window to sum over in seconds (default 0.1s)
        time_offset -- time offset added to each trigger time shift the window (default 0.0)
        max_neuron_id -- the number of neurons
        ''' 
        hist = np.zeros(max_neuron_id) 
        for t_spike in trigger_times:
            ts = t_spike+time_offset
            spikes = self.get_spikes(ts, ts+time_window, max_neuron_id)
            if len(spikes):
                sar = np.array(spikes, dtype=int)[:,1]
                counts = np.bincount(sar, minlength=max_neuron_id)
                hist += counts
        return hist

def isi(spikes):
    '''
    Computes the ISI for spikes 
    '''

    last_spikes = np.zeros(1)
    ISIs = []
    for spike in spikes:
        t,i = spike
        if i>=len(last_spikes):
            last_spikes.resize(i+1)
        if last_spikes[i] > 0:
            ISIs.append(t-last_spikes[i])
        last_spikes[i] = t
    return ISIs

def isi_hist( spikes, *args, **kwargs ):
    pl.hist(isi(spikes), *args, **kwargs )

def cvisi(spikes):
    '''
    Computes the CV ISI for spikes 
    '''

    last_spikes = np.zeros(1)
    sum1 = np.zeros(1)
    sum2 = np.zeros(1)
    nspikes = np.zeros(1)
    for spike in spikes:
        t,i = spike
        if i>=len(last_spikes):
            last_spikes.resize(i+1)
            sum1.resize(i+1)
            sum2.resize(i+1)
            nspikes.resize(i+1)
        if last_spikes[i] > 0:
            isi = t-last_spikes[i]
            sum1[i] += isi
            sum2[i] += isi**2
            nspikes[i] += 1
        last_spikes[i] = t

    cvisi_dist = []
    for i in xrange(len(sum1)):
        if nspikes[i]<2: continue
        mean = sum1[i]/nspikes[i]
        var  = sum2[i]/(nspikes[i]-1)-mean**2
        cvisi_dist.append(np.sqrt(var)/mean)

    return cvisi_dist

def cvisi_hist( spikes, *args, **kwargs ):
    pl.hist(cvisi(spikes), *args, **kwargs )

def spike_counts(spikes):
    '''
    Compute spike count for each neuron
    '''

    nspikes = np.zeros(1)
    for spike in spikes:
        t,i = spike
        if i>=len(nspikes):
            nspikes.resize(i+1)
        nspikes[i] += 1
    return nspikes

def rates(spikes):
    '''
    Compute firing rates for each neuron
    '''

    nspikes = spike_counts(spikes)
    ar = np.array(spikes)
    t_min = ar[:,0].min()
    t_max = ar[:,0].max()
    t_diff = t_max-t_min
    return 1.0*nspikes/t_diff

def rate_hist( spikes, *args, **kwargs ):
    pl.hist(rates(spikes), *args, **kwargs )


def main():
    # running the example program sim_coba_binmon will
    # generate this file
    filenames = ["/tmp/coba.0.e.spk"]
    spkfile = AurynBinarySpikeView(filenames)
    spikes = spkfile.get_spikes()

    rate_hist(spikes, bins=50)
    pl.show()

    # t_start = 0.0
    # t_end   = 0.3
    # n_max = 200

    # spikes = np.array(spkfile.get_spikes(t_start,t_end,max_id=n_max))

    # pl.scatter(spikes[:,0], spikes[:,1])
    # pl.xlabel("Time [s]")
    # pl.ylabel("Neuron ID")
    # pl.show()
    

if __name__ == "__main__":
    main()

