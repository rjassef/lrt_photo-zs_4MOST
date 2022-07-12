import numpy as np

from .runLRT import RunLRT

class GetSEDFits(RunLRT):

    def __init__(self, with_AGN=True, ztype="zspec"):
        super(GetSEDFits,self).__init__(fit_type="SED_fit", with_AGN=with_AGN, ztype=ztype)
        return
