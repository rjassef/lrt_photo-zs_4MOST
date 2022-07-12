import numpy as np

from .runLRT import RunLRT

class GetPhotozs(RunLRT):

    def __init__(self, with_AGN=True):
        super(GetPhotozs,self).__init__(fit_type="zphot", with_AGN=with_AGN, ztype="zphot")
        return
