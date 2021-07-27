import xarray as xr
import numpy as np

""" 
Contains Radiosonde class to get different derived variables (
saturation pressure water vapour, IWV, and specific humidity)
    
Implementation by Andreas Walbroel and Theresa Kiszler, 2021.
"""
class RS_data(object):
    """ Create an instance of a radiosonde data set. 
        Attributes: 
            temp
            press
            rh
            IWV_gg  integrated water vapor following   
            q     specific humidity
    """
    # global constants:
    global M_dv, g
    R_d = 287.04    # gas constant of dry air, in J kg^-1 K^-1
    R_v = 461.5     # gas constant of water vapour, in J kg^-1 K^-1
    M_dv = R_d / R_v # molar mass ratio , in ()
#    e_0 = 611       # saturation water vapour pressure at freezing point (273.15 K), in Pa
#    T0 = 273.15     # freezing temperature, in K
    g = 9.80665     # gravitation acceleration, in m s^-2 (from https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication330e2008.pdf)

    def __init__(self,RS_ncfile):
        """ Initialize radiosonde data
        Parameters:
        -----------
        RS_ncfile : string
            path to netcdf file containing temp, pres and relative humidity
        """
        rs_ds = xr.open_dataset(RS_ncfile)
        self.temp = rs_ds.temp.values
        self.pres = rs_ds.press.values
        self.rh   = rs_ds.rh.values
        self.qi_gg    = np.nan
        self.e_sat_gg = np.nan
        self.IWV_gg  = np.nan

    def get_IWV_gg(self):
        """
        Compute Integrated Water Vapour (also known as precipitable water content)
        out of specific humidity (in kg kg^-1), gravitational constant and air pressure (in Pa).
        Parameters:
        -----------
        q : array of floats
            One dimensional array of specific humidity in kg kg^-1.
        press : array of floats
            One dimensional array of pressure in hPa.
        Returns:
        --------
        IWV : float
            Integrated water vapour kg m-2
        """
        q = self.get_q_gg()
        IWV = self.IWV_gg
        pres = self.pres

        # Check if the Pressure axis is sorted in descending order:
        if np.any((pres[1:] - pres[:-1]) > 0):
            pdb.set_trace()
            raise ValueError("Height axis must be in ascending order to compute the integrated" +
                " water vapour.")
        n_height = len(pres)
        IWV = 0.0
        for k in range(n_height):
            if k == 0:      # bottom of grid
                dp = 0.5*(pres[k+1] - pres[k])      # just a half of a level difference
                IWV = IWV - q[k]*dp
            elif k == n_height-1:   # top of grid
                dp = 0.5*(pres[k] - pres[k-1])      # the other half level difference
                IWV = IWV - q[k]*dp
            else:
                dp = 0.5*(pres[k+1] - pres[k-1])
                IWV = IWV - q[k]*dp
        IWV = IWV / g       # yet had to be divided by gravitational acceleration
        return IWV

    def get_e_sat_gg(self):
        """
        Calculates the saturation pressure over water after Goff and Gratch (1946).
        It is the most accurate that you can get for a temperture range from -90°C to +80°C.
        Source: Smithsonian Tables 1984, after Goff and Gratch 1946
        http://cires.colorado.edu/~voemel/vp.html
        http://hurri.kean.edu/~yoh/calculations/satvap/satvap.html
        e_sat in Pa.
        Parameters:
        -----------
        temp : array of floats
            Array of temperature (in K).
        """
        temp = self.temp
        self.e_sat_gg = 100 * 1013.246 * 10**(-7.90298*(373.16/temp-1) + 5.02808*np.log10(
            373.16/temp) - 1.3816e-7*(10**(11.344*(1-temp/373.16))-1) + 8.1328e-3 * (10**(-3.49149*(373.16/temp-1))-1)) 
        return self.e_sat_gg

    def get_q_gg(self):
        """
        Convert array of relative humidity (between 0 and 1) to specific humidity
        in kg kg^-1.
        Saturation water vapour pressure computation is based on: see e_sat(temp).
        Parameters:
        -----------
        temp : array of floats
            Array of temperature (in K).
        pres : array of floats
            Array of pressure (in Pa).
        relhum : array of floats
            Array of relative humidity (between 0 and 1).
        """
        e_sat_water = self.get_e_sat_gg()
        e = e_sat_water * self.rh
        self.q_gg = M_dv * e / (e*(M_dv - 1) + self.pres)
        return self.q_gg
