# Helpful things

### Radiosonde class
Reading of radiosonde data (tested for 2020 Ny-Alesund data) and directly computing IWV, specific humidity and saturation pressure based on Goff and Gratch 1946.

@Andreas Walbroel, Theresa Kiszler  

Necessary packages: `xarray` and `numpy`. Example usage:  

```
import RS_class

filepath="/path/to/netcdffile.nc"
rs_data = RS_class.RS_data(filepath)

temp = rs_data.temp
IWV = rs_data.get_IWV_gg()

```

### Fortran_cheatsheet
A not yet finished Fortran cheat sheet containing a collection of basic commands which I find useful for my purposes, and possibly for others to.

It is always possible that there may be some mistakes.
