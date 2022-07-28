import numpy as np

from .runLRT import RunLRT

class GetStarFits(RunLRT):

    def __init__(self, stype="MS"):
        super(GetStarFits,self).__init__(fit_type="star_fit", ztype="star", stype=stype)
        return
