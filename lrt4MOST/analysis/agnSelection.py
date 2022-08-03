import numpy as np

class AGNSelection(object):

    def savelist(self, fname, kuse, res, chi2_max, zphot_min_cut):

        kuse = kuse & (res.chi2_agn<chi2_max)
        if zphot_min_cut and self.ztype=="zphot":
            kuse = kuse & (res.z>0.01)

        self.objid_AGN = res.id[kuse]
        ahat_AGN = res.ahat[kuse]
        z_AGN = res.z[kuse]
        ebv_AGN = res.ebv[kuse]
        np.savetxt(fname,np.array([self.objid_AGN, z_AGN, ahat_AGN, ebv_AGN]).T, fmt='%15.0f %10.3f %7.2f %7.2f')

        return

    def is_not_star(self, res, res_stars, select_stars_function):
        kuse = np.ones(res.chi2_agn.shape, dtype=np.bool)
        if res_stars is not None:
            if select_stars_function is None:
                kuse = (kuse) & (res_stars.chi2>res.chi2_agn)
            else:
                kuse = (kuse) & (~select_stars_function(res, res_stars))
        return kuse