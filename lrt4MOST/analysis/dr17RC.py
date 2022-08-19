import numpy as np
from scipy.interpolate import interp1d
import pkg_resources

class DR17RC(object):

    def __init__(self, ttype, data_folder="Cumulative_Arrays"):

        #Save the test type.
        self.ttype = ttype
        self.data_folder = data_folder

        #Read the cumulative number of sources.
        self.stypes = ["QSO", "GALAXY", "STAR"]
        self.bin_edge = dict()
        self.cumsum = dict()
        for stype in self.stypes:

            filename = self.data_folder+"/Cumulative_{}_{}.dat".format(self.ttype, stype)
            stream = pkg_resources.resource_stream(__name__, filename)

            self.rmag_bins = np.array([float(ix) for ix in stream.readline().split()])

            data = np.loadtxt(stream)

            self.bin_edge[stype] = data[:,0]
            if self.ttype=='BIC':
                self.bin_edge[stype] = 10**(self.bin_edge[stype])
            self.cumsum[stype] = data[:,1:]

            stream.close()

        #Set the interpolation functions. 
        self.cumsum_interp = dict()
        for stype in self.stypes:
            self.cumsum_interp[stype] = list()
            for k in range(len(self.rmag_bins)+1):
                self.cumsum_interp[stype].append(interp1d(self.bin_edge[stype], self.cumsum[stype][:,k], fill_value='extrapolate'))

        return

    def C(self, par, rmag):

        #Set the output array
        Cval = np.zeros(len(par))

        #Go per bin of magnitude. 
        for k,rmag_bin in enumerate(self.rmag_bins):
            if k==0:
                rb = -np.inf
                rf = rmag_bin
            else:
                rb = self.rmag_bins[k-1]
                rf = rmag_bin

            #Set the normalization
            norm = np.max(self.cumsum['QSO'][:,k+1])

            #Get the completeness
            cond = (rmag>rb) & (rmag<=rf)
            Cval[cond] = self.cumsum_interp['QSO'][k+1](par[cond])/norm

        #For objects outside of the magnitude ranges we have defined, we will use the entire sample. 
        norm = np.max(self.cumsum['QSO'][:,0])
        cond = (rmag>self.rmag_bins[-1])
        Cval[cond] = self.cumsum_interp['QSO'][0](par[cond])/norm

        return Cval

    def R(self, par, rmag):

        #Set the output array
        Rval = np.zeros(len(par))

        #Go per bin of magnitude. 
        for k,rmag_bin in enumerate(self.rmag_bins):
            if k==0:
                rb = -np.inf
                rf = rmag_bin
            else:
                rb = self.rmag_bins[k-1]
                rf = rmag_bin

            #Get the number of sources in each class. 
            cond = (rmag>rb) & (rmag<=rf)
            NQSO = self.cumsum_interp['QSO'][k+1](par[cond])
            NGAL = self.cumsum_interp['GALAXY'][k+1](par[cond])
            NSTAR = self.cumsum_interp['STAR'][k+1](par[cond])
            Rval[cond] = NQSO/(NQSO+NGAL+NSTAR)

        #For objects outside of the magnitude ranges we have defined, we will use the entire sample. 
        cond = (rmag>self.rmag_bins[-1])
        NQSO = self.cumsum_interp['QSO'][0](par[cond])
        NGAL = self.cumsum_interp['GALAXY'][0](par[cond])
        NSTAR = self.cumsum_interp['STAR'][0](par[cond])
        Rval[cond] = NQSO/(NQSO+NGAL+NSTAR)

        return Rval


