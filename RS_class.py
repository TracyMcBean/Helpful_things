import xarray as xr
import numpy as np

""" 
Contains Radiosonde class to get different derived variables (
saturation pressure water vapour, IWV, and specific humidity)
    
Implementation by Andreas Walbroel and Theresa Kiszler, 2021.

References to saturation water vapour formulae: 
http://cires.colorado.edu/~voemel/vp.html
http://hurri.kean.edu/~yoh/calculations/satvap/satvap.html
"""

class RS_data_hw(object):
    """ Create an instance of a radiosonde data set.
        Computing IWV using the formula by Hyland and Wexler (1983, used by GRUAN)
    """
    # global constants:
    global M_dv, g, R_v
    R_d = 287.04    # gas constant of dry air, in J kg^-1 K^-1
    R_v = 461.5     # gas constant of water vapour, in J kg^-1 K^-1
    M_dv = R_d / R_v # molar mass ratio , in ()
    g = 9.80665     # gravitation acceleration, in m s^-2 (from https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication330e2008.pdf)

    def __init__(self,RS_ncfile):
        """ Initialize radiosonde data
        Parameters:
        -----------
        RS_ncfile : string
            path to netcdf file containing temp, pres and relative humidity
        """
        rs_ds = xr.open_dataset(RS_ncfile)
        self.temp = rs_ds.temp.values     # Temperature in Kelvin
        self.pres = rs_ds.press.values
        self.rh   = rs_ds.rh.values       # relative humidity (0 to 1)
        self.alt  = rs_ds.alt.values      # height in meters
        self.geopot = rs_ds.geopot.values # geopotential height
        self.q_hw    = np.nan    # specific humidity in kg kg^-1
        self.hua_hw = np.nan     # absolute humidity kg m^-3
        self.e_sat_hw = np.nan   # saturation vapour pressure Pa
        self.IWV_hw  = np.nan    # integrated water vapour kg m^-2

    def get_IWV_hw(self):
        """
        Compute Integrated Water Vapour (also known as precipitable water content)
        out of specific humidity (in kg kg^-1), gravitational constant and air pressure (in Pa).

        Returns:
        --------
        IWV : float
            Integrated water vapour kg m-2
        """
        hua = self.get_hua_hw()
        IWV = self.IWV_hw
        geopot = self.geopot
        
        dh = geopot[1:] - geopot[:-1]
        IWV = sum(hua[:-1]*dh)
        
        """ this doesn't work yet 
        prev_geopot = geopot[0]
        for k in range(1,len(geopot)):
            print(geopot[k])
            if geopot[k] < prev_geopot:
                print("skipping value where geopot is lower than previous one.")
            else:
                dh = geopot[k] - prev_geopot
                IWV += hua[k]*dh
                prev_geopot = geopot[k-1]  # So that we only use values where the height is increasing
        """

        return IWV

    def get_e_sat_hw(self):
        """
        Compute saturation water vapour pressure following Hyland and Wexler (1983)
        """
        temp = self.temp
        self.e_sat_hw = np.exp(  -0.58002206e4/temp + 0.13914993e1  - 0.48640239e-1*temp + 0.41764768e-4*temp**2 - 0.14452093e-7 *temp**3 + 0.65459673e1 * np.log(temp))
        return self.e_sat_hw

    def get_q_hw(self):
        """
        Convert array of relative humidity (between 0 and 1) to specific humidity
        in kg kg^-1.
        Saturation vapour pressure calculated by Hyland and Wexler (1983)
        """
        
        e_sat_water = self.get_e_sat_hw()
        e = e_sat_water * self.rh
        self.q_hw = M_dv * e / (e*(M_dv - 1) + self.pres)
        return self.q_hw

    def get_hua_hw(self):
        """
        Get absolute humidity using saturation vapour pressure calculated by Hyland and Wexler (1983).
        """
        
        e_env = self.rh * self.get_e_sat_hw()
        self.hua_hw = e_env / (R_v*self.temp)
        return self.hua_hw
    
    
class RS_data_gg(object):
    """ 
    Create an instance of a radiosonde data set. 
    Computing IWV using the formula by Goff and Gratch (1946)
    """
    # global constants:
    global M_dv, g, R_v
    R_d = 287.04    # gas constant of dry air, in J kg^-1 K^-1
    R_v = 461.5     # gas constant of water vapour, in J kg^-1 K^-1
    M_dv = R_d / R_v # molar mass ratio , in ()
    g = 9.80665     # gravitation acceleration, in m s^-2 (from https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication330e2008.pdf)

    def __init__(self,RS_ncfile):
        """ Initialize radiosonde data
        Parameters:
        -----------
        RS_ncfile : string
            path to netcdf file containing temp, pres and relative humidity
        """
        rs_ds = xr.open_dataset(RS_ncfile)
        self.temp = rs_ds.temp.values     # Temperature in Kelvin
        self.pres = rs_ds.press.values
        self.rh   = rs_ds.rh.values       # relative humidity (0 to 1)
        self.alt  = rs_ds.alt.values      # height in meters
        self.geopot = rs_ds.geopot.values # geopotential height
        self.q_gg    = np.nan    # specific humidity in kg kg^-1
        self.hua_gg = np.nan     # absolute humidity kg m^-3
        self.e_sat_gg = np.nan   # saturation vapour pressure Pa
        self.IWV_gg  = np.nan    # integrated water vapour kg m^-2

    def get_IWV_gg(self):
        """
        Compute Integrated Water Vapour (also known as precipitable water content)
        out of specific humidity (in kg kg^-1), gravitational constant and air pressure (in Pa).
        Returns
        --------
        IWV : float
            Integrated water vapour kg m-2
        """
        
        hua = self.get_hua_gg()
        IWV = self.IWV_gg
        geopot = self.geopot
        
        dh = geopot[1:] - geopot[:-1]
        IWV = sum(hua[:-1]*dh)
        
        return IWV

    def get_e_sat_gg(self):
        """
        Calculates the saturation pressure over water after Goff and Gratch (1946).
        It is the most accurate that you can get for a temperture range from -90°C to +80°C.
        Source: Smithsonian Tables 1984, after Goff and Gratch 1946
        http://cires.colorado.edu/~voemel/vp.html
        http://hurri.kean.edu/~yoh/calculations/satvap/satvap.html
        e_sat in Pa.
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
        """
        e_sat_water = self.get_e_sat_gg()
        e = e_sat_water * self.rh
        self.q_gg = M_dv * e / (e*(M_dv - 1) + self.pres)
        return self.q_gg
    
    def get_hua_gg(self):
        """
        Get absolute humidity using saturation vapour pressure calculated by Hyland and Wexler (1983).
        """
        
        e_env = self.rh * self.get_e_sat_gg()
        self.hua_gg = e_env / (R_v*self.temp)
        return self.hua_gg
    
