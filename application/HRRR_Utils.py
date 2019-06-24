
from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from netCDF4 import num2date
import numpy as np
from metpy.units import units
import scipy.ndimage as ndimage
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
# =============================================================================
# RETRIEVE RAP AND HRRR DATA
# =============================================================================
RAP = 'http://thredds-jetstream.unidata.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_20km/latest.xml'
HRRR= 'http://thredds-jetstream.unidata.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/latest.xml'
GFS = 'http://thredds-jetstream.unidata.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p5deg/latest.xml'

DATA = TDSCatalog(GFS)
BEST_DATA = list(DATA.datasets.values())[0]
NCSS_DATA = NCSS(BEST_DATA.access_urls['NetcdfSubset'])

NOW = datetime.utcnow()
LATEST_DATA = NCSS_DATA.query().time(NOW).accept('netcdf4')
# =============================================================================
# UPPER-AIR VARIABLES
# =============================================================================
# 250: JET STREAM, GEOPOTENTIAL HEIGHT, POTENTIAL VORTICITY, IRROTATIONAL WIND
def 250hPa_GFS_jet_stream_SLP(lon_west, lon_east, lat_south, lat_north):
       
def 250hPa_GFS_jet_stream_jet_dyn(lon_west, lon_east, lat_south, lat_north):
# 500: VORTICITY, GEOPOTENTIAL HEIGHT, VORTICITY ADVECTION
def 500hPa_GFS_vorticity(lon_west, lon_east, lat_south, lat_north):
       LATEST_DATA.variables('Geopotential_height_isobaric', 'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric').add_lonlat()
       LATEST_DATA.lonlat_box(lon_west, lon_east, lat_south, lat_north)
       LATEST_DATA.vertical_level(50000)
       DATA = NCSS_DATA.get_data(LATEST_DATA)
       DTIME = DATA.variables['Geopotential_height_isobaric'].dimensions[0]
       DLAT = DATA.variables['Geopotential_height_isobaric'].dimensions[2]
       DLON = DATA.variables['Geopotential_height_isobaric'].dimensions[3]
       LAT = DATA.variables[DLAT][:]
       LON = DATA.variables[DLON][:]
       TIMES = DATA.variables[DTIME]
       VTIMES = num2date(TIMES[:], TIMES.units)
       H500 = DATA.variables['Geopotential_height_isobaric'][0, 0, :, :]
       U500 = DATA.variables['u-component_of_wind_isobaric'][0, 0, :, :]*units('m/s')
       V500 = DATA.variables['v-component_of_wind_isobaric'][0, 0, :, :]*units('m/s')
       DTIME = DATA.variables['Geopotential_height_isobaric'].dimensions[:]
       f = mpcalc.coriolis_parameter(np.deg2rad(LAT)).to(units('1/sec'))
       DX, DY = mpcalc.lat_lon_grid_spacing(LON, LAT)
       VORT500 = mpcalc.vorticity(U500, V500, DX, DY, dim_order='YX')
       VORT500 = (VORT500*(units('1/s'))) 


# 700: Q-VECTORS+CONVERGENCE, GEOPOTENTIAL HEIGHT, POTENTIAL TEMPERATURE
# 850: GEOPOTENTIAL HEIGHT, WINDS, TEMP-ADVECTION, FRONTOGENESIS
# 925: GEOPOTENTIAL HEIGHT, WINDS, TEMP-ADVECTION, FRONTOGENESIS
# SFC: MSLP, WIND, TEMPERATURE
# =============================================================================
# THERMODYNAMIC VARIABLES
# =============================================================================
#MUCAPE: Pre-computed in HRRR
#LATEST_DATA.variables('Convective_available_potential_energy_surface').add_lonlat()
#LATEST_DATA.lonlat_box(280, 295, 35, 50)
#DATA = NCSS_DATA.get_data(LATEST_DATA)
#SCAPE = DATA.variables['Convective_available_potential_energy_surface'][:]
#LAT = DATA.variables['lat'][:]
#LON = DATA.variables['lon'][:]
#MLCAPE: Pre-computed in HRRR
#SBCAPE: Pre-computed in HRRR
#LCL: Pre-computed in HRRR
#LI: Pre-computed in HRRR
#EL: Pre-computed in HRRR  
#LAPSE RATES: Needs Temps various levels to compute. Easy.
#CIN: Pre-computed in HRRR
# =============================================================================
# WIND VARIABLES
# =============================================================================
# SFC-1km SHEAR: Pre-computed in HRRR
# SFC-6km SHEAR: Pre-computed in HRRR
# SFC-1km HELICITY: Pre-computed in HRRR
# SFC-3km HELICITY: Pre-computed in HRRR
# UPSLOPE FLOW: Requires 80-m Wind and Terrain Height[requires surface pressure and temps]
# =============================================================================
# COMPOSITE INDICES
# =============================================================================
# SUPERCELL
# FIXED-LAYER SIG TORNADO
# LARGE HAIL
# CRAVEN BROOKS
# =============================================================================
# PRECIPITATION
# =============================================================================
# SIMULATED REFLECTIVITY
# PRECIPITATION ACCUMULATION
# PRECIPITABLE WATER
# MOISTURE FLUX/MOISTURE FLUX CONVERGENCE
# SNOWFALL

datacrs = ccrs.PlateCarree()
plotcrs = ccrs.Miller(central_longitude=-95.0)

# Make a grid of lat/lon values to use for plotting with Basemap.
lons, lats = np.meshgrid(np.squeeze(LON), np.squeeze(LAT))

fig = plt.figure(1, figsize=(12., 13.))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02],
                       bottom=.07, top=.99, hspace=0.01, wspace=0.01)

ax = plt.subplot(gs[0], projection=plotcrs)
ax.set_title('500-hPa Geopotential Heights (m)', loc='left')
ax.set_title('VALID: {}'.format(VTIMES[0]), loc='right')

#   ax.set_extent([west long, east long, south lat, north lat])
ax.set_extent([230, 300, 25, 55], ccrs.PlateCarree())
ax.coastlines('50m', edgecolor='black', linewidth=0.5)
ax.add_feature(cfeature.STATES, linewidth=0.5)

clev250 = np.arange(4800, 6600, 60)
cs = ax.contour(lons, lats, H500, clev250, colors='k', linewidths=1.0, linestyles='solid', transform=datacrs)
plt.clabel(cs, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
           rightside_up=True, use_clabeltext=True)

clevsped250 = np.arange(0, 55, 5)
cmap = plt.cm.get_cmap('YlOrBr')
cf = ax.contourf(lons, lats, np.squeeze(VORT500*10**5), clevsped250, cmap=cmap, transform=datacrs)
cax = plt.subplot(gs[1])
cbar = plt.colorbar(cf, cax=cax, orientation='horizontal', extend='max', extendrect=True)

gs.tight_layout(fig)
plt.show()  

  