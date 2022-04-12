#load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
#load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
#load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
#load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
#***********************************************
#begin
#***********************************************

#import numpy,sys,os
#import Nio
#import Ngl

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import matplotlib
matplotlib.style.use('ggplot')
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.font_manager import FontProperties
import numpy as np
import os
from wrf import getvar, to_np, ll_to_xy, vertcross, smooth2d, CoordPair, get_basemap, latlon_coords, interplevel
import sys
from datetime import datetime, timedelta

def model(directory, inputfilename):

  n = 1
  T2_newa = np.zeros(0)
  Q2_newa = np.zeros(0)
  T2_newb = np.zeros(0)
  Q2_newb = np.zeros(0)
  psfc_newa = np.zeros(0)
  psfc_newb = np.zeros(0)
  Q_newa = np.zeros(0)
  Q_newb = np.zeros(0)
  TH2_newa = np.zeros(0)
  TH2_newb = np.zeros(0)
  T2_newc = np.zeros(0)
  Q2_newc = np.zeros(0) 
  T2_newd = np.zeros(0)
  Q2_newd = np.zeros(0) 
  psfc_newc = np.zeros(0)
  psfc_newd = np.zeros(0)
  Q_newc = np.zeros(0)
  Q_newd = np.zeros(0)
  TH2_newc = np.zeros(0)
  TH2_newd = np.zeros(0)
  inputFileList = open(inputfilename, "r")
  for filename in inputFileList:
    if not filename.rstrip().endswith(".py"):
      if not filename.rstrip().endswith(".ncl"):
        nci = nc.Dataset(directory+filename.rstrip())
        coorda = ll_to_xy(nci, 50.8210, -115.2141, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #BNS
        coordb = ll_to_xy(nci, 50.8260, -115.1983, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #POW 
        coordd = ll_to_xy(nci, 50.9568, -115.2044, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #FIS 
        coordc = ll_to_xy(nci, 50.9441, -115.1389, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #Hay 
        locXa = coorda[0]
        locYa = coorda[1]
        locXb = coordb[0]
        locYb = coordb[1]
        locXc = coordc[0]
        locYc = coordc[1]
        locXd = coordd[0]
        locYd = coordd[1]       
       
        T2a = nci['T2'][0,locYa,locXa]
        T2_newa = np.append(T2_newa,T2a)
        q2a = nci['Q2'][0,locYa,locXa]
        Q2_newa = np.append(Q2_newa,q2a)
        T2b = nci['T2'][0,locYb,locXb]
        T2_newb = np.append(T2_newb,T2b)
        q2b = nci['Q2'][0,locYb,locXb]
        Q2_newb = np.append(Q2_newb,q2b)
        psfca = nci['PSFC'][0,locYa,locXa]
        psfc_newa = np.append(psfc_newa,psfca)
        psfcb = nci['PSFC'][0,locYb,locXb]
        psfc_newb = np.append(psfc_newb,psfcb)
        qa = nci['QVAPOR'][0,2,locYa,locXa]
        Q_newa = np.append(Q_newa,qa)
        qb = nci['QVAPOR'][0,2,locYb,locXb]
        Q_newb = np.append(Q_newb,qb)
        TH2a = nci['TH2'][0,locYa,locXa]
        TH2_newa = np.append(TH2_newa,TH2a)
        TH2b = nci['TH2'][0,locYb,locXb]
        TH2_newb = np.append(TH2_newb,TH2b)
        
        T2c = nci['T2'][0,locYc,locXc]
        T2_newc = np.append(T2_newc,T2c)
        q2c = nci['Q2'][0,locYc,locXc]
        Q2_newc = np.append(Q2_newc,q2c)
        T2d = nci['T2'][0,locYd,locXd]
        T2_newd = np.append(T2_newd,T2d)
        q2d = nci['Q2'][0,locYd,locXd]
        Q2_newd = np.append(Q2_newd,q2d)
        psfcc = nci['PSFC'][0,locYc,locXc]
        psfc_newc = np.append(psfc_newc,psfcc)
        psfcd = nci['PSFC'][0,locYd,locXd]
        psfc_newd = np.append(psfc_newd,psfcd)
        qc = nci['QVAPOR'][0,2,locYc,locXc]
        Q_newc = np.append(Q_newc,qc)
        qd = nci['QVAPOR'][0,2,locYd,locXd]
        Q_newd = np.append(Q_newd,qd)
        TH2c = nci['TH2'][0,locYc,locXc]
        TH2_newc = np.append(TH2_newc,TH2c)
        TH2d = nci['TH2'][0,locYd,locXd]
        TH2_newd = np.append(TH2_newd,TH2d)

        n = n + 1 
  return T2_newa, Q2_newa, T2_newb, Q2_newb, psfc_newa, psfc_newb, Q_newa, Q_newb, TH2_newa, TH2_newb, T2_newc, Q2_newc, T2_newd, Q2_newd, psfc_newc, psfc_newd, Q_newc, Q_newd, TH2_newc, TH2_newd 
#
data_obs_b = pd.read_csv('/glade/work/mina/data_python_sodar/POWER_2016.csv')
data_obs_a = pd.read_csv('/glade/work/mina/data_python_sodar/BNS_2016.csv')
data_obs_c = pd.read_csv('/glade/work/mina/data_python_sodar/Marmot_hay_2016.csv')
data_obs_d = pd.read_csv('/glade/work/mina/data_python_sodar/Marmot_fis_2016.csv')

inputFileName = "inputdata"
directory = "/WRFV3.7.1_ERA_PBL/WRFV3/run/output/domain4/"
T2_new_pbla, Q2_new_pbla, T2_new_pblb, Q2_new_pblb, psfc_newa, psfc_newb, Q_newa, Q_newb, TH2_newa, TH2_newb, T2_new_pblc, Q2_new_pblc, T2_new_pbld, Q2_new_pbld, psfc_newc, psfc_newd, Q_newc, Q_newd, TH2_newc, TH2_newd = model(directory, inputFileName)
directory = "/WRFV3.7.1_ERA_LES/WRFV3/run/output/domain4/"
T2_new_lesa, Q2_new_lesa, T2_new_lesb, Q2_new_lesb, psfc_newa, psfc_newb, Q_newa, Q_newb, TH2_newa, TH2_newb, T2_new_lesc, Q2_new_lesc, T2_new_lesd, Q2_new_lesd, psfc_newc, psfc_newd, Q_newc, Q_newd, TH2_newc, TH2_newd = model(directory, inputFileName)


T2_new_pbla = T2_new_pbla - 273.15
T2_new_lesa = T2_new_lesa - 273.15
T2_new_pblb = T2_new_pblb - 273.15
T2_new_lesb = T2_new_lesb - 273.15
T2_new_pblc = T2_new_pblc - 273.15
T2_new_lesc = T2_new_lesc - 273.15
T2_new_pbld = T2_new_pbld - 273.15
T2_new_lesd = T2_new_lesd - 273.15


nci = nc.Dataset("wrfout_d04_2016-07-18_16:45:00")

coorda = ll_to_xy(nci, 50.8210, -115.2141, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #BON 3.5 or 5m  used 5m
coordb = ll_to_xy(nci, 50.8260, -115.1983, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #POW 2 or 5m    used 5m
coordd = ll_to_xy(nci, 50.9568, -115.2044, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #FIS 2.6m
coordc = ll_to_xy(nci, 50.9441, -115.1389, timeidx=0, squeeze=True, meta=True, stagger=None, as_int=True)  #Hay 7m
locXa = coorda[0]
locYa = coorda[1]
locXb = coordb[0]
locYb = coordb[1]


data_t_a = pd.read_csv('BNS.d04.T.jul18.TS',sep="\s+", header = None,low_memory=False)
data_t_b = pd.read_csv('POW.d04.T.jul18.TS',sep="\s+", header = None,low_memory=False)
data_t_c = pd.read_csv('HAY.d04.T.jul18.TS',sep="\s+", header = None,low_memory=False)
data_t_d = pd.read_csv('FIS.d04.T.jul18.TS',sep="\s+", header = None,low_memory=False)
time = pd.read_csv('POW.d04.time.jul18',sep="\s+", header = None,low_memory=False)
data_t_a_g = pd.read_csv('BNS.d04.T.jul18.TS.LESGF',sep="\s+", header = None,low_memory=False)
data_t_b_g = pd.read_csv('POW.d04.T.jul18.TS.LESGF',sep="\s+", header = None,low_memory=False)
data_t_c_g = pd.read_csv('HAY.d04.T.jul18.TS.LESGF',sep="\s+", header = None,low_memory=False)
data_t_d_g = pd.read_csv('FIS.d04.T.jul18.TS.LESGF',sep="\s+", header = None,low_memory=False)
data_t_a_s = pd.read_csv('BNS.d04.T.jul18.TS.shade',sep="\s+", header = None,low_memory=False)
data_t_b_s = pd.read_csv('POW.d04.T.jul18.TS.shade',sep="\s+", header = None,low_memory=False)
data_t_c_s = pd.read_csv('HAY.d04.T.jul18.TS.shade',sep="\s+", header = None,low_memory=False)
data_t_d_s = pd.read_csv('FIS.d04.T.jul18.TS.shade',sep="\s+", header = None,low_memory=False)
time_s = pd.read_csv('BNS.d04.time.jul18.shade',sep="\s+", header = None,low_memory=False)

chunk_size_ave = 21684     # data pounts for 15 min LES averaging

groups = [data_t_a[x:x+chunk_size_ave] for x in range(0, len(data_t_a), chunk_size_ave)]
Ta_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Ta_ave[i] = np.nanmean(group)

groups = [data_t_b[x:x+chunk_size_ave] for x in range(0, len(data_t_b), chunk_size_ave)]
Tb_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Tb_ave[i] = np.nanmean(group)

groups = [data_t_c[x:x+chunk_size_ave] for x in range(0, len(data_t_c), chunk_size_ave)]
Tc_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Tc_ave[i] = np.nanmean(group)

groups = [data_t_d[x:x+chunk_size_ave] for x in range(0, len(data_t_d), chunk_size_ave)]
Td_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Td_ave[i] = np.nanmean(group)

groups_time = [time[x:x+chunk_size_ave] for x in range(0, len(time), chunk_size_ave)]
time_ave= np.zeros(len(groups_time))
for i in range(len(groups_time)):
  group_time = np.array(groups_time[i])
  time_ave[i] = np.nanmean(group_time)

groups = [data_t_a_s[x:x+chunk_size_ave] for x in range(0, len(data_t_a_s), chunk_size_ave)]
Ta_s_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Ta_s_ave[i] = np.nanmean(group)

groups = [data_t_b_s[x:x+chunk_size_ave] for x in range(0, len(data_t_b_s), chunk_size_ave)]
Tb_s_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Tb_s_ave[i] = np.nanmean(group)

groups = [data_t_c_s[x:x+chunk_size_ave] for x in range(0, len(data_t_c_s), chunk_size_ave)]
Tc_s_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Tc_s_ave[i] = np.nanmean(group)

groups = [data_t_d_s[x:x+chunk_size_ave] for x in range(0, len(data_t_d_s), chunk_size_ave)]
Td_s_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Td_s_ave[i] = np.nanmean(group)

groups = [data_t_a_g[x:x+chunk_size_ave] for x in range(0, len(data_t_a_g), chunk_size_ave)]
Ta_g_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Ta_g_ave[i] = np.nanmean(group)

groups = [data_t_b_g[x:x+chunk_size_ave] for x in range(0, len(data_t_b_g), chunk_size_ave)]
Tb_g_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Tb_g_ave[i] = np.nanmean(group)

groups = [data_t_c_g[x:x+chunk_size_ave] for x in range(0, len(data_t_c_g), chunk_size_ave)]
Tc_g_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Tc_g_ave[i] = np.nanmean(group)

groups = [data_t_d_g[x:x+chunk_size_ave] for x in range(0, len(data_t_d_g), chunk_size_ave)]
Td_g_ave= np.zeros(len(groups))
for i in range(len(groups)):
  group = np.array(groups[i])
  Td_g_ave[i] = np.nanmean(group)

groups_time = [time_s[x:x+chunk_size_ave] for x in range(0, len(time_s), chunk_size_ave)]
time_s_ave= np.zeros(len(groups_time))
for i in range(len(groups_time)):
  group_time = np.array(groups_time[i])
  time_s_ave[i] = np.nanmean(group_time)

time = [5,5.16,5.33,5.5,5.66,5.83,6,6.16,6.33,6.5,6.66,6.83,7,7.16,7.33,7.5,7.66,7.83,8,8.16,8.33,8.5,8.66,8.83,9,9.16,9.33,9.5,9.66,9.83,10,10.16,10.33,10.5,10.66,10.83,11,11.16,11.33,11.5,11.66,11.83,12,12.16,12.33,12.5,12.66,12.83,13,13.16,13.33,13.5,13.66,13.83,14,14.16,14.33,14.5,14.66,14.83,15,15.16,15.33,15.5,15.66,15.83,16,16.16,16.33,16.5,16.66,16.83,17,17.16,17.33,17.5,17.66,17.83,18,18.16,18.33,18.5,18.66,18.83,19,19.16,19.33,19.5,19.66,19.83,20,20.16,20.33,20.5,20.66,20.83,21]

data_time = data_obs_a['time'][1:67].to_numpy()
data_time = np.mean(data_time.reshape(-1,2),axis=1)

Ta_ave = Ta_ave - 273.15
Tb_ave = Tb_ave - 273.15
Tc_ave = Tc_ave - 273.15
Td_ave = Td_ave - 273.15
Ta_s_ave = Ta_s_ave - 273.15
Tb_s_ave = Tb_s_ave - 273.15
Tc_s_ave = Tc_s_ave - 273.15
Td_s_ave = Td_s_ave - 273.15
Ta_g_ave = Ta_g_ave - 273.15
Tb_g_ave = Tb_g_ave - 273.15
Tc_g_ave = Tc_g_ave - 273.15
Td_g_ave = Td_g_ave - 273.15


fontP = FontProperties()
fontP.set_size('small')
plt.figure()
plt.rcParams.update(plt.rcParamsDefault)
plt.subplot(2,2,1)
axes = plt.gca()
axes.set_xlim([5,21])
axes.set_ylim([0,25])
plt.ylabel('Air temperature ($^\circ$C)')
plt.plot(data_obs_a['time'][1:65],data_obs_a['AirTemp_Avg'][1:65], 'k-o', markersize=4, markeredgecolor='k', linewidth = 1, label="OBS")
plt.plot(time_ave-1.0,Ta_ave, '-', markersize=4, color = 'green', linewidth = 1.2, label="LESLF")
plt.plot(time_s_ave-7.0,Ta_s_ave, '-', markersize=4, color = 'lime', linewidth = 1.2, label="LESLF_shade")
plt.plot(time,T2_new_pbla[0:97],':', markersize=4, color = 'maroon', linewidth = 1.2, label="PBL")
plt.plot(time_ave[0:72]-1.0,Ta_g_ave, '--', markersize=4, color = 'darkblue', linewidth = 1.2, label="LESGF")
plt.text(20,22,'(a)',weight="bold")
plt.legend(loc=2, borderaxespad=0.,fancybox=True,shadow=False,fontsize='small')

axes = plt.gca()
plt.subplot(2,2,2)
axes = plt.gca()
axes.set_xlim([5,21])
axes.set_ylim([0,25])
plt.ylabel('Air temperature ($^\circ$C)')
plt.plot(data_obs_a['time'][1:65],data_obs_b['AirTemp_Avg'][1:65], 'k-o', markersize=4, markeredgecolor='k', linewidth = 1, label="OBS")
plt.plot(time_ave-1.0,Tb_ave, '-', markersize=4, color = 'green', linewidth = 1.2, label="LESLF")
plt.plot(time_s_ave-7.0,Tb_s_ave, '-', markersize=4, color = 'lime', linewidth = 1.2, label="LESLF_shade")
plt.plot(time,T2_new_pblb[0:97],':', markersize=4, color = 'maroon', linewidth = 1.2, label="PBL")
plt.plot(time_ave[0:72]-1.0,Tb_g_ave, '--', markersize=4, color = 'darkblue', linewidth = 1.2, label="LESGF")
plt.text(20,22,'(b)',weight="bold")
plt.legend(loc=2, borderaxespad=0.,fancybox=True,shadow=False,fontsize='small')


axes = plt.gca()
plt.subplot(2,2,3)
axes = plt.gca()
axes.set_xlim([5,21])
axes.set_ylim([0,25])
plt.ylabel('Air temperature ($^\circ$C)')
plt.xlabel('Local Time (h)')
plt.plot(data_obs_a['time'][1:65],data_obs_c['AirTemp_Avg'][1:65], 'k-o', markersize=4, markeredgecolor='k', linewidth = 1, label="OBS")
plt.plot(time_ave-1.0,Tc_ave, '-', markersize=4, color = 'green', linewidth = 1.2, label="LESLF")
plt.plot(time_s_ave-7.0,Tc_s_ave, '-', markersize=4, color = 'lime', linewidth = 1.2, label="LESLF_shade")
plt.plot(time_ave[0:72]-1.0,Tc_g_ave, '--', markersize=4, color = 'darkblue', linewidth = 1.2, label="LESGF")
plt.plot(time,T2_new_pblc[0:97],':', markersize=4, color = 'maroon', linewidth = 1.2, label="PBL")
plt.text(20,22,'(c)',weight="bold")

axes = plt.gca()
plt.subplot(2,2,4)
axes = plt.gca()
axes.set_xlim([5,21])
axes.set_ylim([0,25])
plt.ylabel('Air temperature ($^\circ$C)')
plt.xlabel('Local Time (h)')
plt.plot(data_obs_a['time'][1:65],data_obs_d['AirTemp_Avg'][1:65], 'k-o', markersize=4, markeredgecolor='k', linewidth = 1, label="OBS")
plt.plot(time_ave-1.0,Td_ave, '-', markersize=4, color = 'green', linewidth = 1.2, label="LESLF")
plt.plot(time_s_ave-7.0,Td_s_ave, '-', markersize=4, color = 'lime', linewidth = 1.2, label="LESLF_shade")
plt.plot(time_ave[0:72]-1.0,Td_g_ave, '--', markersize=4, color = 'darkblue', linewidth = 1.2, label="LESGF")
plt.plot(time,T2_new_pbld[0:97],':', markersize=4, color = 'maroon', linewidth = 1.2, label="PBL")
plt.text(20,22,'(d)',weight="bold")
plt.savefig('BNS_POW_HAY_FIS_20160718_Temp.png')   # save the figure to file
plt.show()
plt.close()

