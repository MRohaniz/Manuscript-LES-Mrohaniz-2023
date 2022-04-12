#***********************************************
#***********************************************

import numpy as np
import scipy 
from scipy import signal
from scipy.fftpack import fftfreq
import pandas as pd
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
matplotlib.style.use('ggplot')
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.font_manager import FontProperties
from matplotlib import rc
import netCDF4 as nc
import os
from wrf import getvar, to_np, ll_to_xy, vertcross, smooth2d, CoordPair, get_basemap, latlon_coords, interplevel
import sys
from datetime import datetime, timedelta

def model(directory, inputfilename):

  n = 1
  speed_newa = np.zeros(0)
  direc_newa = np.zeros(0)
  speed_newb = np.zeros(0)
  direc_newb = np.zeros(0)
  inputFileList = open(inputfilename, "r")
  for filename in inputFileList:
    if not filename.rstrip().endswith(".py"):
      if not filename.rstrip().endswith(".ncl"):
        nci = nc.Dataset(directory+filename.rstrip())
        coorda = ll_to_xy(nci, 50.8210, -115.2141, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #location of BNS 
        coordb = ll_to_xy(nci, 50.8260, -115.1983, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #location of POW 
        locXa = coorda[0]
        locYa = coorda[1]
        locXb = coordb[0]
        locYb = coordb[1]
        Uaa = nci['U10'][0,locYa,locXa]
        Vaa = nci['V10'][0,locYa,locXa]
        Ubb = nci['U10'][0,locYb,locXb]
        Vbb = nci['V10'][0,locYb,locXb]

        speed10a = (Uaa**2+Vaa**2)**0.5
        USTa = nci['UST'][0,locYa,locXa]
        speeda = speed10a - ((USTa/0.4) * np.log(10/5))

        speed10b = (Ubb**2+Vbb**2)**0.5
        USTb = nci['UST'][0,locYb,locXb]
        speedb = speed10b - ((USTb/0.4) * np.log(10/5))
        speed_newa = np.append(speed_newa,speeda)
        speed_newb = np.append(speed_newb,speedb)

        cosalphaa = nci['COSALPHA'][0,locYa,locXa]
        sinalphaa = nci['SINALPHA'][0,locYa,locXa]
        Uaa = Uaa*cosalphaa - Vaa*sinalphaa
        Vaa = Vaa*cosalphaa + Uaa*sinalphaa
        r2d = 45.0/np.arctan(1.0)
        direca = np.arctan2(Uaa, Vaa) * r2d + 180
        direc_newa = np.append(direc_newa,direca)

        cosalphab = nci['COSALPHA'][0,locYb,locXb]
        sinalphab = nci['SINALPHA'][0,locYb,locXb]
        Ubb = Ubb*cosalphab - Vbb*sinalphab
        Vbb = Vbb*cosalphab + Ubb*sinalphab
        direcb = np.arctan2(Ubb, Vbb) * r2d + 180
        direc_newb = np.append(direc_newb,direcb)
        n = n + 1

  return speed_newa, direc_newa,speed_newb, direc_newb

inputFileName = "inputdata"
directory = "WRFV3.7.1_ERA_PBL/WRFV3/run/output/domain4/"
speed_new_pbla, direc_new_pbla,speed_new_pblb, direc_new_pblb = model(directory, inputFileName)
directory = "WRFV3.7.1_ERA_LES/WRFV3/run/output/domain4/"
speed_new_lesa, direc_new_lesa, speed_new_lesb, direc_new_lesb = model(directory, inputFileName)


nci = nc.Dataset("wrfout_d04_2016-07-18_16:45:00")

coorda = ll_to_xy(nci, 50.8210, -115.2141, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #BNS 
coordb = ll_to_xy(nci, 50.8260, -115.1983, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #POW 
locXa = coorda[0]
locYa = coorda[1]
locXb = coordb[0]
locYb = coordb[1]

data_u_a = pd.read_csv('BNS.d04.U.jul18.TS',sep="\s+", header = None,low_memory=False)
data_v_a = pd.read_csv('BNS.d04.V.jul18.TS',sep="\s+", header = None,low_memory=False)
data_u_b = pd.read_csv('POW.d04.U.jul18.TS',sep="\s+", header = None,low_memory=False)
data_v_b = pd.read_csv('POW.d04.V.jul18.TS',sep="\s+", header = None,low_memory=False)
time = pd.read_csv('POW.d04.time.jul18',sep="\s+", header = None,low_memory=False)
data_u_a_g = pd.read_csv('BNS.d04.U.jul18.TS.LESGF',sep="\s+", header = None,low_memory=False)
data_v_a_g = pd.read_csv('BNS.d04.V.jul18.TS.LESGF',sep="\s+", header = None,low_memory=False)
data_u_b_g = pd.read_csv('POW.d04.U.jul18.TS.LESGF',sep="\s+", header = None,low_memory=False)
data_v_b_g = pd.read_csv('POW.d04.V.jul18.TS.LESGF',sep="\s+", header = None,low_memory=False)
data_u_a_s = pd.read_csv('BNS.d04.U.jul18.TS.shade',sep="\s+", header = None,low_memory=False)
data_v_a_s = pd.read_csv('BNS.d04.V.jul18.TS.shade',sep="\s+", header = None,low_memory=False)
data_u_b_s = pd.read_csv('POW.d04.U.jul18.TS.shade',sep="\s+", header = None,low_memory=False)
data_v_b_s = pd.read_csv('POW.d04.V.jul18.TS.shade',sep="\s+", header = None,low_memory=False)
time_s = pd.read_csv('BNS.d04.time.jul18.shade',sep="\s+", header = None,low_memory=False)


Ua = data_u_a
Ub = data_u_b
Va = data_v_a
Vb = data_v_b
Ua_g = data_u_a_g
Ub_g = data_u_b_g
Va_g = data_v_a_g
Vb_g = data_v_b_g
Ua_s = data_u_a_s
Ub_s = data_u_b_s
Va_s = data_v_a_s
Vb_s = data_v_b_s

speed10a = (Ua**2+Va**2)**0.5
USTa = nci['UST'][0,locYa,locXa]
speeda = np.abs(speed10a - ((USTa/0.4) * np.log(10/5)))

speed10b = (Ub**2+Vb**2)**0.5
USTb = nci['UST'][0,locYb,locXb]
speedb = np.abs(speed10b - ((USTb/0.4) * np.log(10/5)))

speed10ag = (Ua_g**2+Va_g**2)**0.5
USTa = nci['UST'][0,locYa,locXa]
speeda_g = np.abs(speed10ag - ((USTa/0.4) * np.log(10/5)))

speed10bg = (Ub_g**2+Vb_g**2)**0.5
USTb = nci['UST'][0,locYb,locXb]
speedb_g = np.abs(speed10bg - ((USTb/0.4) * np.log(10/5)))

speeda_s = (Ua_s**2+Va_s**2)**0.5
USTa = nci['UST'][0,locYa,locXa]

speed10bs = (Ub_s**2+Vb_s**2)**0.5
speedb_s = (Ub_s**2+Vb_s**2)**0.5
USTb = nci['UST'][0,locYb,locXb]

cosalphaa = nci['COSALPHA'][0,locYa,locXa]
sinalphaa = nci['SINALPHA'][0,locYa,locXa]
cosalphab = nci['COSALPHA'][0,locYb,locXb]
sinalphab = nci['SINALPHA'][0,locYb,locXb]

r2d = 45.0/np.arctan(1.0)
direca = np.arctan2(Ua, Va) * r2d + 180
direcb = np.arctan2(Ub, Vb) * r2d + 180
direca_g = np.arctan2(Ua_g, Va_g) * r2d + 180
direcb_g = np.arctan2(Ub_g, Vb_g) * r2d + 180
direca_s = np.arctan2(Ua_s, Va_s) * r2d + 180
direcb_s = np.arctan2(Ub_s, Vb_s) * r2d + 180


chunk_size_ave = 21684     # number of data points for 15 min averaging for LES timeseries data

groups = [speeda[x:x+chunk_size_ave] for x in range(0, len(speeda), chunk_size_ave)]
speeda_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  speeda_ave[i] = np.nanmean(group)

groups = [speedb[x:x+chunk_size_ave] for x in range(0, len(speedb), chunk_size_ave)]
speedb_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  speedb_ave[i] = np.nanmean(group)


groups = [direca[x:x+chunk_size_ave] for x in range(0, len(direca), chunk_size_ave)]
direca_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  direca_ave[i] = np.nanmean(group)

groups = [direcb[x:x+chunk_size_ave] for x in range(0, len(direcb), chunk_size_ave)]
direcb_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  direcb_ave[i] = np.nanmean(group)

groups_time = [time[x:x+chunk_size_ave] for x in range(0, len(time), chunk_size_ave)]
time_ave= np.zeros(len(groups_time))
for i in range(len(groups_time)):
  group_time = np.array(groups_time[i])
  time_ave[i] = np.nanmean(group_time)


groups = [speeda_s[x:x+chunk_size_ave] for x in range(0, len(speeda_s), chunk_size_ave)]
speeda_ave_s= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  speeda_ave_s[i] = np.nanmean(group)

groups = [speedb_s[x:x+chunk_size_ave] for x in range(0, len(speedb_s), chunk_size_ave)]
speedb_ave_s= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  speedb_ave_s[i] = np.nanmean(group)


groups = [direca_s[x:x+chunk_size_ave] for x in range(0, len(direca_s), chunk_size_ave)]
direca_ave_s= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  direca_ave_s[i] = np.nanmean(group)

groups = [direcb_s[x:x+chunk_size_ave] for x in range(0, len(direcb_s), chunk_size_ave)]
direcb_ave_s= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  direcb_ave_s[i] = np.nanmean(group)

groups_time = [time_s[x:x+chunk_size_ave] for x in range(0, len(time_s), chunk_size_ave)]
time_ave_s= np.zeros(len(groups_time))
for i in range(len(groups_time)):
  group_time = np.array(groups_time[i])
  time_ave_s[i] = np.nanmean(group_time)

groups = [speeda_g[x:x+chunk_size_ave] for x in range(0, len(speeda_g), chunk_size_ave)]
speeda_ave_g= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  speeda_ave_g[i] = np.nanmean(group)

groups = [speedb_g[x:x+chunk_size_ave] for x in range(0, len(speedb_g), chunk_size_ave)]
speedb_ave_g= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  speedb_ave_g[i] = np.nanmean(group)


groups = [direca_g[x:x+chunk_size_ave] for x in range(0, len(direca_g), chunk_size_ave)]
direca_ave_g= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  direca_ave_g[i] = np.nanmean(group)

groups = [direcb_g[x:x+chunk_size_ave] for x in range(0, len(direcb_g), chunk_size_ave)]
direcb_ave_g= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  direcb_ave_g[i] = np.nanmean(group)


data_obs_a = pd.read_csv('/glade/work/mina/data_python_sodar/BNS_2016.csv')
data_obs_b = pd.read_csv('/glade/work/mina/data_python_sodar/POWER_2016.csv')

time = [5,5.16,5.33,5.5,5.66,5.83,6,6.16,6.33,6.5,6.66,6.83,7,7.16,7.33,7.5,7.66,7.83,8,8.16,8.33,8.5,8.66,8.83,9,9.16,9.33,9.5,9.66,9.83,10,10.16,10.33,10.5,10.66,10.83,11,11.16,11.33,11.5,11.66,11.83,12,12.16,12.33,12.5,12.66,12.83,13,13.16,13.33,13.5,13.66,13.83,14,14.16,14.33,14.5,14.66,14.83,15,15.16,15.33,15.5,15.66,15.83,16,16.16,16.33,16.5,16.66,16.83,17,17.16,17.33,17.5,17.66,17.83,18,18.16,18.33,18.5,18.66,18.83,19,19.16,19.33,19.5,19.66,19.83,20,20.16,20.33,20.5,20.66,20.83,21]


fontP = FontProperties()
fontP.set_size('small')
plt.figure()
plt.rcParams.update(plt.rcParamsDefault)
plt.subplot(2,2,1)
axes = plt.gca()
axes.set_xlim([5,21])
axes.set_ylim([0,10])
axes = plt.gca()
plt.ylabel('WS (m/s)')
plt.plot(data_obs_a['time'][0:65],data_obs_a['WindSpeed_S_WVT'][0:65], 'k-o', markersize=4, markeredgecolor='k', linewidth = 1, label="OBS")
plt.plot(time_ave-1.0,speeda_ave, '-', markersize=4, color = 'green', linewidth = 1.2, label="LESLF")
plt.plot(time_ave_s-7.0,speeda_ave_s, '-', markersize=4, color = 'lime', linewidth = 1.2, label="LESLF_shade")
plt.plot(time_ave[0:72]-1.0,speeda_ave_g, '--', markersize=4, color = 'darkblue', linewidth = 1.2, label="LESGF")
plt.plot(time,speed_new_pbla[0:97], ':', markersize=4, color = 'maroon', linewidth = 1.2, label="PBL")
plt.text(20,9,'(a)',weight="bold")
plt.legend(loc=2, borderaxespad=0.,fancybox=True,shadow=False,fontsize='small')


axes = plt.gca()
plt.subplot(2,2,2)
axes = plt.gca()
axes.set_xlim([5,21])
axes.set_ylim([0,10])
plt.ylabel('WS (m/s)')
plt.plot(data_obs_a['time'][0:65],data_obs_b['WindSpeed_S_WVT'][0:65], 'k-o', markersize=4, markeredgecolor='k', linewidth = 1, label="OBS")
plt.plot(time_ave-1.0,speedb_ave, '-', markersize=4, color = 'green', linewidth = 1.2, label="LESLF")
plt.plot(time_ave_s-7.0,speedb_ave_s, '-', markersize=4, color = 'lime', linewidth = 1.2, label="LESLF_shade")
plt.plot(time_ave[0:72]-1.0,speedb_ave_g, '--', markersize=4, color = 'darkblue', linewidth = 1.2, label="LESGF")
plt.plot(time,speed_new_pblb[0:97], ':', markersize=4, color = 'maroon', linewidth = 1.2, label="PBL")
plt.text(20,9,'(b)',weight="bold")
plt.legend(loc=2, borderaxespad=0.,fancybox=True,shadow=False,fontsize='small')


axes = plt.gca()
plt.subplot(2,2,3)
axes = plt.gca()
axes.set_xlim([5,21])
axes.set_ylim(-50,360)
plt.ylabel('WD')
plt.xlabel('Local Time (h)')
plt.plot(data_obs_a['time'][0:65],data_obs_a['WindDir_D1_WVT'][0:65], 'k-o', markersize=4, markeredgecolor='k', linewidth = 1, label="OBS")  
plt.plot(time_ave-1.0,direca_ave, '-', markersize=4, color = 'green', linewidth = 1.2, label="LESLF")
plt.plot(time_ave_s-7.0,direca_ave_s, '-', markersize=4, color = 'lime', linewidth = 1.2, label="LESLF_shade")
plt.plot(time_ave[0:72]-1.0,direca_ave_g, '--', markersize=4, color = 'darkblue', linewidth = 1.2, label="LESGF")
plt.plot(time,direc_new_pbla[0:97], ':', markersize=4, color = 'maroon', linewidth = 1.2, label="PBL")

axes = plt.gca()
plt.subplot(2,2,4)
axes = plt.gca()
axes.set_xlim([5,21])
axes.set_ylim(-50,360)
plt.ylabel('WD')
plt.xlabel('Local Time (h)')
plt.plot(data_obs_a['time'][0:65],data_obs_b['WindDir_D1_WVT'][0:65], 'k-o', markersize=4, markeredgecolor='k', linewidth = 1, label="OBS")
plt.plot(time_ave-1.0,direcb_ave, '-', markersize=4, color = 'green', linewidth = 1.2, label="LESLF")
plt.plot(time_ave_s-7.0,direcb_ave_s, '-', markersize=4, color = 'lime', linewidth = 1.2, label="LESLF_shade")
plt.plot(time_ave[0:72]-1.0,direcb_ave_g, '--', markersize=4, color = 'darkblue', linewidth = 1.2, label="LESGF")
plt.plot(time,direc_new_pblb[0:97], ':', markersize=4, color = 'maroon', linewidth = 1.2, label="PBL")
plt.savefig('Figure_BNS_POW_wind_shade.png')   # save the figure to file
plt.show()
plt.close()

sys.exit(0)



