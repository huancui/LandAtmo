#!/usr/bin/env python
# coding: utf-8

# ### Script to count the number of events in each 5 degree boxes in order to check how it matches with other studies

# In[ ]:


import os
import sys
import numpy as np
from netCDF4 import Dataset
import csv
import string
from scipy import stats

year1 = 1997
year2 = 2019

months        = ['04','05','06','07','08']
days          = [30,31,30,31,31]
monthstr      = ['apr','may','jun','jul','aug']
july_pentads  = [1,6,11,16,21,26,31]

path1     = '/pic/projects/flood/CONUS_simulation/'
path5     = '/pic/projects/flood/CONUS_simulation/LDASOUT_LandAtmo/SMgradient/'

mask_data = np.genfromtxt(path1 + 'CentralEastUS_mask.csv', delimiter=',')

filename1 = path1 + 'geo_em.d01.nc'
data1     = Dataset(filename1, 'r', format='NETCDF4')
lat1      = data1.variables['CLAT']
lon1      = data1.variables['CLONG']
landmask  = data1.variables['LANDMASK']     # 0->water,  1->land
height    = data1.variables['HGT_M']
print(lon1.shape)
lon11     = lon1[0,:,:]
lat11     = lat1[0,:,:]
land11    = landmask[0,:,:]
height11   = height[0,:,:]
print(lon11.shape)


# In[ ]:


#-----define ~5 degree blocks-----
box_size         = 6
block_size       = 20

land_1         = np.mean(np.reshape(land11, (    lon11.shape[0],                     int(lon11.shape[1]/box_size), box_size)), axis=2)
land_2         = np.mean(np.reshape(land_1, (int(lon11.shape[0]/box_size), box_size, int(lon11.shape[1]/box_size)          )), axis=1)

lon_1          = np.mean(np.reshape(lon11, (    lon11.shape[0],                     int(lon11.shape[1]/box_size), box_size)), axis=2)
lon_2          = np.mean(np.reshape(lon_1, (int(lon11.shape[0]/box_size), box_size, int(lon11.shape[1]/box_size)          )), axis=1)
print(lon_2.shape)

lat_1           = np.mean(np.reshape(lat11, (    lon11.shape[0],                     int(lon11.shape[1]/box_size), box_size)), axis=2)
lat_2           = np.mean(np.reshape(lat_1, (int(lon11.shape[0]/box_size), box_size, int(lon11.shape[1]/box_size)          )), axis=1)

#-----in order to reduce the number of precip samples, add 2 more criterion-----
#height >1500
height_1        = np.mean(np.reshape(height11, (    lon11.shape[0],                     int(lon11.shape[1]/box_size), box_size)), axis=2)
height_2        = np.mean(np.reshape(height_1, (int(lon11.shape[0]/box_size), box_size, int(lon11.shape[1]/box_size)          )), axis=1)


centralUS_mask  = np.zeros((int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))

lat_low    = 25.
lat_up     = 50.
lon_left   = -120.
lon_right  = -80.

for jdim in range(centralUS_mask.shape[0]):
    for idim in range(centralUS_mask.shape[1]):
        if (lat_2[jdim, idim] >= lat_low and lat_2[jdim, idim] <= lat_up):
            if (lon_2[jdim, idim] >= lon_left and lon_2[jdim, idim] <= lon_right):
                centralUS_mask[jdim,idim] = 1.0

#----define 5 degree boxes----
box_degree  = 5
box_degree2 = 1
lat_blocks_start       = np.zeros((int((lat_up-lat_low)/box_degree2), int((lon_right-lon_left)/box_degree2)))
lon_blocks_start       = np.zeros((int((lat_up-lat_low)/box_degree2), int((lon_right-lon_left)/box_degree2)))
lat_blocks_end         = np.zeros((int((lat_up-lat_low)/box_degree2), int((lon_right-lon_left)/box_degree2)))
lon_blocks_end         = np.zeros((int((lat_up-lat_low)/box_degree2), int((lon_right-lon_left)/box_degree2)))
centralUS_blocks       = np.zeros((int((lat_up-lat_low)/box_degree2), int((lon_right-lon_left)/box_degree2)))
lat_blocks             = np.zeros((int((lat_up-lat_low)/box_degree2), int((lon_right-lon_left)/box_degree2)))
lat_blocks_center      = np.zeros((int((lat_up-lat_low)/box_degree2), int((lon_right-lon_left)/box_degree2)))
lon_blocks_center      = np.zeros((int((lat_up-lat_low)/box_degree2), int((lon_right-lon_left)/box_degree2)))
print(lat_blocks.shape)

for jdim in range(lat_blocks_start.shape[0]):
    for idim in range(lat_blocks_start.shape[1]):
        lat_blocks_center[jdim, idim]= lat_low  + jdim*box_degree2 + box_degree2/2.0
        lon_blocks_center[jdim, idim]= lon_left + idim*box_degree2 + box_degree2/2.0
        lat_blocks_start[jdim, idim] = lat_blocks_center[jdim, idim] - box_degree/2.0
        lat_blocks_end[jdim, idim]   = lat_blocks_center[jdim, idim] + box_degree/2.0
        lon_blocks_start[jdim, idim] = lon_blocks_center[jdim, idim] - box_degree/2.0
        lon_blocks_end[jdim, idim]   = lon_blocks_center[jdim, idim] + box_degree/2.0
        


# In[ ]:


#----define functions----
def define_block(lat00, lon00):
    print(centralUS_mask.shape)
    print(lat_2.shape)
    map1 = np.zeros(centralUS_mask.shape)
    lat01 = lat00 + box_degree
    lon01 = lon00 + box_degree
    map1  = np.where(np.logical_and(lat_2>=lat00, lat_2<lat01), 1, 0)
    map1  = np.where(np.logical_and(lon_2>=lon00, lon_2<lon01), map1, 0)
    return map1

import datetime
def get_datenum(y1,m1,d1,m2):
    time1 = datetime.datetime(y1,m1,d1)
    time0 = datetime.datetime(y1,m2+4, 1)
    return int((time1-time0).days + (y1-1997)*days[m2])

def calculate_bootstrap(data_event, data_control):
    sample_num              = len(data_event)
    data_all                = np.concatenate([data_event, data_control])
    nonevent_bootstrap      = np.zeros((1000, sample_num))
    nonevent_bootstrap_mean = np.zeros(1000)
    for iboot in range(1000):
        bootstrapi                     = np.random.choice(data_all, size=sample_num)
        rest_data                      = np.delete(data_all, bootstrapi)
        nonevent_bootstrap[iboot,:]    = bootstrapi
        nonevent_bootstrap_mean[iboot] = np.mean(bootstrapi)-np.mean(rest_data)
    return stats.percentileofscore(nonevent_bootstrap_mean, (np.mean(data_event)-np.mean(data_control)))
    


# In[ ]:


imonth                     = 3
max_event                  = 3000
max_nonevent               = days[imonth]*(year2-year1) 

PP_allyears       = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
PPsur_allyears    = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
PPamount_allyears = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
MCSmask_allyears  = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
AMmask_allyears   = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
year_allyears     = np.zeros(days[imonth]*(year2-year1))

#-----define matrix to store percentile data-----
SM1m_allyears                      = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
SM1m_allyears_prime                = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
SM1m_delta_PERCENTILES             = np.zeros((3, lat_blocks.shape[0], lat_blocks.shape[1]))
SM1m_anomaly_PERCENTILES           = np.zeros((3, lat_blocks.shape[0], lat_blocks.shape[1]))
SM1mTR_allyears                    = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
SM1mTR_allyears_prime              = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
SM1mTR_delta_PERCENTILES           = np.zeros((3, lat_blocks.shape[0], lat_blocks.shape[1]))
SM1mTR_anomaly_PERCENTILES         = np.zeros((3, lat_blocks.shape[0], lat_blocks.shape[1]))
SM1mNTR_allyears                   = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
SM1mNTR_allyears_prime             = np.zeros((days[imonth]*(year2-year1), int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
SM1mNTR_delta_PERCENTILES          = np.zeros((3, lat_blocks.shape[0], lat_blocks.shape[1]))
SM1mNTR_anomaly_PERCENTILES        = np.zeros((3, lat_blocks.shape[0], lat_blocks.shape[1]))
SM1m_delta_PERCENTILES[:,:,:]      = np.nan
SM1m_anomaly_PERCENTILES[:,:,:]    = np.nan
SM1mTR_delta_PERCENTILES[:,:,:]    = np.nan
SM1mTR_anomaly_PERCENTILES[:,:,:]  = np.nan
SM1mNTR_delta_PERCENTILES[:,:,:]   = np.nan
SM1mNTR_anomaly_PERCENTILES[:,:,:] = np.nan

for iyear in range(year1, year2):
    #-----top 1m soil moisture-----
    filename1    = path5 + 'Obtain_9amSM1m_'+str(iyear)+'_24km120km.nc'
    data1        = Dataset(filename1,'r',format='NETCDF4')
    soilm1m       = data1.variables['soilm1']
    soilm1mtr     = data1.variables['soilm1_mcs']
    soilm1mntr    = data1.variables['soilm1_nonmcs']
    
    ind1         = (iyear-year1)*days[imonth]
    ind2         = (iyear-year1+1)*days[imonth]
    ind3         = sum(days[0:imonth])
    ind4         = ind3 + days[imonth]
    print('iyear,  ind1, ind2: ', iyear,  ind1, ind2)
    print('imonth, ind3, ind4: ', imonth, ind3, ind4)
    SM1m_allyears[ind1:ind2,:,:]       = soilm1m[ind3:ind4,:,:]
    SM1mTR_allyears[ind1:ind2,:,:]     = soilm1mtr[ind3:ind4,:,:]
    SM1mNTR_allyears[ind1:ind2,:,:]    = soilm1mntr[ind3:ind4,:,:]
    
    filename2    = path5 + 'Obtain_AfternoonPP_maxind_' + str(iyear) + '_24km120km_iterate_filterRX1.nc'
    data2        = Dataset(filename2,'r',format='NETCDF4')
    max_pp       = data2.variables['MAXPP_ind']
    pp_amount    = data2.variables['PP_amount']
    mask_morning = data2.variables['MAX_morningmask']  # 1->no mornong rain, 2->morning rain
    mask_mcs     = data2.variables['MAX_mcsmask']      # 1->non-MCS, 2->MCS
    max_ppsur    = data2.variables['MAXPP_sur'] 
    PP_allyears[ind1:ind2,:,:]       = max_pp[ind3:ind4,:,:]
    PPamount_allyears[ind1:ind2,:,:] = pp_amount[ind3:ind4,:,:]
    MCSmask_allyears[ind1:ind2,:,:]  = mask_mcs[ind3:ind4,:,:]
    AMmask_allyears[ind1:ind2,:,:]   = mask_morning[ind3:ind4,:,:]
    PPsur_allyears[ind1:ind2,:,:]    = max_ppsur[ind3:ind4,:,:]
    year_allyears[ind1:ind2]         = iyear
print('end loading all years for month:', imonth)

#-----remove seasonal cycle-----
for iday in range(days[imonth]):
    SM1m_mean          = np.nanmean(SM1m_allyears[iday:(year2-year1)*days[imonth]:days[imonth],:,:],axis=0)
    SM1mTR_mean        = np.nanmean(SM1mTR_allyears[iday:(year2-year1)*days[imonth]:days[imonth],:,:],axis=0)
    SM1mNTR_mean       = np.nanmean(SM1mNTR_allyears[iday:(year2-year1)*days[imonth]:days[imonth],:,:],axis=0)
    for iyear in range(year1, year2):
        time_ind = (iyear-year1)*days[imonth]+iday
        SM1m_allyears_prime[time_ind,:,:]       = np.subtract(SM1m_allyears[time_ind,:,:], SM1m_mean)
        SM1mTR_allyears_prime[time_ind,:,:]     = np.subtract(SM1mTR_allyears[time_ind,:,:], SM1mTR_mean)
        SM1mNTR_allyears_prime[time_ind,:,:]    = np.subtract(SM1mNTR_allyears[time_ind,:,:], SM1mNTR_mean)
print('time_ind at the last step after removing seasonal cycle:', time_ind)

#-----calculate SM preference for each 5-degree box-----
for jdim in range(lat_blocks.shape[0]):
    for idim in range(lat_blocks.shape[1]):
        pixel_5degree = define_block(lat_blocks_start[jdim, idim],lon_blocks_start[jdim, idim])
        
        #-----create matrices to store event data-----
        block_event_data         = np.zeros((max_event, 8))
        block_nonevent_data      = np.zeros((max_event, max_nonevent))
        block_event_data[:,:]    = np.nan
        block_nonevent_data[:,:] = np.nan
        event_ind                = -1
        
        EVENT_SM1m_delta         = np.zeros(max_event)
        EVENT_SM1m_max           = np.zeros(max_event)
        EVENT_MCSmask            = np.zeros(max_event)
        NONEVENT_SM1m_delta      = np.zeros((max_event, max_nonevent))
        NONEVENT_SM1m_max        = np.zeros((max_event, max_nonevent))
        EVENT_SM1mTR_delta          = np.zeros(max_event)
        EVENT_SM1mTR_max            = np.zeros(max_event)
        EVENT_SM1mNTR_delta         = np.zeros(max_event)
        EVENT_SM1mNTR_max           = np.zeros(max_event)
        NONEVENT_SM1mTR_delta       = np.zeros((max_event, max_nonevent))
        NONEVENT_SM1mTR_max         = np.zeros((max_event, max_nonevent))
        NONEVENT_SM1mNTR_delta      = np.zeros((max_event, max_nonevent))
        NONEVENT_SM1mNTR_max        = np.zeros((max_event, max_nonevent))
        
        EVENT_SM1m_delta[:]         = np.nan
        EVENT_SM1m_max[:]           = np.nan
        EVENT_MCSmask[:]            = np.nan
        NONEVENT_SM1m_delta[:,:]    = np.nan
        NONEVENT_SM1m_max[:,:]      = np.nan
        EVENT_SM1mTR_delta[:]       = np.nan
        EVENT_SM1mTR_max[:]         = np.nan
        NONEVENT_SM1mTR_delta[:,:]  = np.nan
        NONEVENT_SM1mTR_max[:,:]    = np.nan
        EVENT_SM1mNTR_delta[:]      = np.nan
        EVENT_SM1mNTR_max[:]        = np.nan
        NONEVENT_SM1mNTR_delta[:,:] = np.nan
        NONEVENT_SM1mNTR_max[:,:]   = np.nan
        
        
        for iyear in range(year1, year2):
            filename1     = path5 + 'Afternoon_PP_max_latind_'+ str(iyear) + months[imonth] + '_24km120km_iterate_filterRX1.txt'
            filename2     = path5 + 'Afternoon_PP_max_lonind_'+ str(iyear) + months[imonth] + '_24km120km_iterate_filterRX1.txt'
            filename3     = path5 + 'Afternoon_PP_min_latind_'+ str(iyear) + months[imonth] + '_24km120km_iterate_filterRX1.txt'
            filename4     = path5 + 'Afternoon_PP_min_lonind_'+ str(iyear) + months[imonth] + '_24km120km_iterate_filterRX1.txt'
            filename5     = path5 + 'Afternoon_PP_datetime_'  + str(iyear) + months[imonth] + '_24km120km_iterate_filterRX1.txt'
            
            f1            = open(filename5, 'r')
            datetime_list = []
            for line in f1.read().split('\n'):
                datetime_list.append(line)
            f2              = open(filename1, 'r')
            max_latind_list = []
            for line in f2.read().split('\n'):
                max_latind_list.append(float(line))
            f3              = open(filename2, 'r')
            max_lonind_list = []
            for line in f3.read().split('\n'):
                max_lonind_list.append(float(line))
            f4              = open(filename3, 'r')
            min_latind_list = []
            for line in f4.readlines():
                line1       = line.strip()
                list        = []
                for num in line1.split(','):
                    list.append(float(num))
                min_latind_list.append(list)
            f5              = open(filename4, 'r')
            min_lonind_list = []
            for line in f5.readlines():
                line1       = line.strip()
                list        = []
                for num in line1.split(','):
                    list.append(float(num))
                min_lonind_list.append(list)
            #----finish reading files----
            
            #==============================================
            #========going through each event==============
            for ii in range(len(datetime_list)):      
                itime    = datetime_list[ii]
                ind_day  = get_datenum(int(itime[0:4]), int(itime[4:6]),int(itime[6:8]), imonth)
                imax_lat = int(max_latind_list[ii])
                imax_lon = int(max_lonind_list[ii])
                imin_lat = min_latind_list[ii]
                imin_lon = min_lonind_list[ii]
                #----tell if within the current 5by5 degree box or not----
                if (pixel_5degree[imax_lat, imax_lon] > 0.5):
                    #-----read in morning mask and MCS mask later-----
                    imax_am     = AMmask_allyears[ind_day, int(imax_lat), int(imax_lon)]
                    imax_mcs    = MCSmask_allyears[ind_day, int(imax_lat), int(imax_lon)]
                    imax_id     = PP_allyears[ind_day,int(imax_lat), int(imax_lon)]
                    imax_times  = PP_allyears[:,int(imax_lat), int(imax_lon)]
                    
                    if (np.isnan(imax_id) ):        # grid not associated with max rain
                        print('ERROR! grid not associated with max rain')
                    else:
                        event_ind                       = event_ind + 1
                        sm1m_max                        = SM1m_allyears_prime[:, int(imax_lat), int(imax_lon)]
                        sm1mtr_max                      = SM1mTR_allyears_prime[:, int(imax_lat), int(imax_lon)]
                        sm1mntr_max                     = SM1mNTR_allyears_prime[:, int(imax_lat), int(imax_lon)]
                        sm1m_mins                       = np.zeros((SM1m_allyears.shape[0], len(imin_lat)))
                        sm1mtr_mins                     = np.zeros((SM1m_allyears.shape[0], len(imin_lat)))
                        sm1mntr_mins                    = np.zeros((SM1m_allyears.shape[0], len(imin_lat)))
                        for jj in range(len(imin_lat)):
                            sm1m_mins[:,jj]             = SM1m_allyears_prime[:, int(imin_lat[jj]), int(imin_lon[jj])]
                            sm1mtr_mins[:,jj]           = SM1mTR_allyears_prime[:, int(imin_lat[jj]), int(imin_lon[jj])]
                            sm1mntr_mins[:,jj]          = SM1mNTR_allyears_prime[:, int(imin_lat[jj]), int(imin_lon[jj])]
                        sm1m_min_ave                    = np.mean(sm1m_mins, axis=1)
                        sm1m_delta                      = np.subtract(sm1m_max, sm1m_min_ave)
                        sm1mtr_min_ave                  = np.mean(sm1mtr_mins, axis=1)
                        sm1mtr_delta                    = np.subtract(sm1mtr_max, sm1mtr_min_ave)
                        sm1mntr_min_ave                 = np.mean(sm1mntr_mins, axis=1)
                        sm1mntr_delta                   = np.subtract(sm1mntr_max, sm1mntr_min_ave)
                        
                        #-------------------------------------------------
                        #   tell which are non-event years
                        imax_year_mask        = np.zeros(days[imonth]*(year2-year1))
                        nonevent_year_mask    = np.zeros(days[imonth]*(year2-year1))
                        nonevent_year_mask[:] = np.nan
                        for jj in range(len(imax_times)):
                            if (imax_times[jj] >= 0):
                                year_event   = year_allyears[jj]
                                ind1         = int((year_event-year1)*days[imonth])
                                ind2         = int((year_event-year1+1)*days[imonth])
                                imax_year_mask[ind1:ind2] =  imax_year_mask[ind1:ind2] +1.0  # year range associated with event
                                
                        for jj in range(len(imax_times)):
                            if (imax_year_mask[jj] == 0.0): #non-event years
                                imax_sur_am = np.where(PPsur_allyears[ind_day, :,:] == imax_id, AMmask_allyears[jj, :, :], np.nan)
                                if (np.nanmax(imax_sur_am) > 2.5):    # 2->morning rain; 1->no morning rain
                                    nonevent_year_mask[jj] = 0.0
                                else:
                                    nonevent_year_mask[jj] = 1.0
                        print('event samples: ',      np.count_nonzero(imax_year_mask>0.5))
                        print('non-event samples: ',  np.count_nonzero(nonevent_year_mask>0.5))
                        multievent_flag = 0
                        if (max(imax_year_mask) > 1.5):
                            print('multiple events in the same year and month!')
                            for kk in np.unique(year_allyears[imax_year_mask>1.5]):
                                if (kk == iyear):
                                    multievent_flag = 2 
                        
                        #----save the information of each event-----
                        block_event_data[event_ind, 0] = int(itime[0:4])
                        block_event_data[event_ind, 1] = int(itime[4:6])
                        block_event_data[event_ind, 2] = int(itime[6:8])
                        block_event_data[event_ind, 3] = imax_id
                        block_event_data[event_ind, 4] = imax_lat
                        block_event_data[event_ind, 5] = imax_lon
                        block_event_data[event_ind, 6] = multievent_flag
                        block_event_data[event_ind, 7] = np.count_nonzero(nonevent_year_mask>0.5)
                        EVENT_MCSmask[event_ind]       = imax_mcs
                        EVENT_SM1m_delta[event_ind]    = sm1m_delta[ind_day]
                        EVENT_SM1m_max[event_ind]      = sm1m_max[ind_day]
                        NONEVENT_SM1m_delta[event_ind,:]= np.where(nonevent_year_mask > 0.5, sm1m_delta, np.nan)
                        NONEVENT_SM1m_max[event_ind,:]  = np.where(nonevent_year_mask > 0.5, sm1m_max, np.nan)
                        EVENT_SM1mTR_delta[event_ind]       = sm1mtr_delta[ind_day]
                        EVENT_SM1mTR_max[event_ind]         = sm1mtr_max[ind_day]
                        NONEVENT_SM1mTR_delta[event_ind,:]  = np.where(nonevent_year_mask > 0.5, sm1mtr_delta, np.nan)
                        NONEVENT_SM1mTR_max[event_ind,:]    = np.where(nonevent_year_mask > 0.5, sm1mtr_max, np.nan)
                        EVENT_SM1mNTR_delta[event_ind]      = sm1mntr_delta[ind_day]
                        EVENT_SM1mNTR_max[event_ind]        = sm1mntr_max[ind_day]
                        NONEVENT_SM1mNTR_delta[event_ind,:] = np.where(nonevent_year_mask > 0.5, sm1mntr_delta, np.nan)
                        NONEVENT_SM1mNTR_max[event_ind,:]   = np.where(nonevent_year_mask > 0.5, sm1mntr_max, np.nan)
                        
        #-----calculate the percentiles for each block-----
        #(1)exclude cases with multievent_flag > 1.5
        if np.count_nonzero(~np.isnan(EVENT_SM1m_delta)) > 25:
            for ievent in range(np.count_nonzero(~np.isnan(EVENT_SM1m_delta))):
                if (block_event_data[ievent, 6] > 1.5): # this year has multiple event for the same grid
                    EVENT_SM1m_delta[ievent]      = np.nan
                    EVENT_SM1m_max[ievent]        = np.nan
                    NONEVENT_SM1m_delta[ievent,:] = np.nan
                    NONEVENT_SM1m_max[ievent,:]   = np.nan
                    EVENT_SM1mTR_delta[ievent]       = np.nan
                    EVENT_SM1mTR_max[ievent]         = np.nan
                    NONEVENT_SM1mTR_delta[ievent,:]  = np.nan
                    NONEVENT_SM1mTR_max[ievent,:]    = np.nan
                    EVENT_SM1mNTR_delta[ievent]      = np.nan
                    EVENT_SM1mNTR_max[ievent]        = np.nan
                    NONEVENT_SM1mNTR_delta[ievent,:] = np.nan
                    NONEVENT_SM1mNTR_max[ievent,:]   = np.nan
            print('event_ind, ievent:', event_ind, ievent)

            event_datax        = EVENT_SM1m_delta[~np.isnan(EVENT_SM1m_delta)]
            nonevent_datax     = NONEVENT_SM1m_delta[~np.isnan(NONEVENT_SM1m_delta)]
            event_datay        = EVENT_SM1m_max[~np.isnan(EVENT_SM1m_max)]
            nonevent_datay     = NONEVENT_SM1m_max[~np.isnan(NONEVENT_SM1m_max)]
            if (len(event_datax) > 25):
                percentile_x   = calculate_bootstrap(event_datax, nonevent_datax)
                SM1m_delta_PERCENTILES[0, jdim, idim] = percentile_x
                percentile_x   = calculate_bootstrap(event_datay, nonevent_datay)
                SM1m_anomaly_PERCENTILES[0, jdim, idim] = percentile_x
            
            event_datax        = EVENT_SM1mTR_delta[~np.isnan(EVENT_SM1mTR_delta)]
            nonevent_datax     = NONEVENT_SM1mTR_delta[~np.isnan(NONEVENT_SM1mTR_delta)]
            event_datay        = EVENT_SM1mTR_max[~np.isnan(EVENT_SM1mTR_max)]
            nonevent_datay     = NONEVENT_SM1mTR_max[~np.isnan(NONEVENT_SM1mTR_max)]
            if (len(event_datax) > 25):
                percentile_x   = calculate_bootstrap(event_datax, nonevent_datax)
                SM1mTR_delta_PERCENTILES[0, jdim, idim] = percentile_x
                percentile_x   = calculate_bootstrap(event_datay, nonevent_datay)
                SM1mTR_anomaly_PERCENTILES[0, jdim, idim] = percentile_x
                
            event_datax        = EVENT_SM1mNTR_delta[~np.isnan(EVENT_SM1mNTR_delta)]
            nonevent_datax     = NONEVENT_SM1mNTR_delta[~np.isnan(NONEVENT_SM1mNTR_delta)]
            event_datay        = EVENT_SM1mNTR_max[~np.isnan(EVENT_SM1mNTR_max)]
            nonevent_datay     = NONEVENT_SM1mNTR_max[~np.isnan(NONEVENT_SM1mNTR_max)]
            if (len(event_datax) > 25):
                percentile_x   = calculate_bootstrap(event_datax, nonevent_datax)
                SM1mNTR_delta_PERCENTILES[0, jdim, idim] = percentile_x
                percentile_x   = calculate_bootstrap(event_datay, nonevent_datay)
                SM1mNTR_anomaly_PERCENTILES[0, jdim, idim] = percentile_x

            #(2) exclude mcs or non-mcs events
            EVENT_SM1m_delta_copy       = np.copy(EVENT_SM1m_delta)
            EVENT_SM1m_max_copy         = np.copy(EVENT_SM1m_max)
            NONEVENT_SM1m_delta_copy    = np.copy(NONEVENT_SM1m_delta)
            NONEVENT_SM1m_max_copy      = np.copy(NONEVENT_SM1m_max)
            EVENT_SM1mTR_delta_copy     = np.copy(EVENT_SM1mTR_delta)
            EVENT_SM1mTR_max_copy       = np.copy(EVENT_SM1mTR_max)
            NONEVENT_SM1mTR_delta_copy  = np.copy(NONEVENT_SM1mTR_delta)
            NONEVENT_SM1mTR_max_copy    = np.copy(NONEVENT_SM1mTR_max)
            EVENT_SM1mNTR_delta_copy    = np.copy(EVENT_SM1mNTR_delta)
            EVENT_SM1mNTR_max_copy      = np.copy(EVENT_SM1mNTR_max)
            NONEVENT_SM1mNTR_delta_copy = np.copy(NONEVENT_SM1mNTR_delta)
            NONEVENT_SM1mNTR_max_copy   = np.copy(NONEVENT_SM1mNTR_max)
            #----exclude non-mcs events first----
            for ievent in range(event_ind+1):
                if (EVENT_MCSmask[ievent] < 1.5):  # an MCS event
                    EVENT_SM1m_delta[ievent]         = np.nan
                    EVENT_SM1m_max[ievent]           = np.nan
                    NONEVENT_SM1m_delta[ievent,:]    = np.nan
                    NONEVENT_SM1m_max[ievent,:]      = np.nan
                    EVENT_SM1mTR_delta[ievent]       = np.nan
                    EVENT_SM1mTR_max[ievent]         = np.nan
                    NONEVENT_SM1mTR_delta[ievent,:]  = np.nan
                    NONEVENT_SM1mTR_max[ievent,:]    = np.nan
                    EVENT_SM1mNTR_delta[ievent]      = np.nan
                    EVENT_SM1mNTR_max[ievent]        = np.nan
                    NONEVENT_SM1mNTR_delta[ievent,:] = np.nan
                    NONEVENT_SM1mNTR_max[ievent,:]   = np.nan
            print('event_ind, ievent:', event_ind, ievent)
            event_datax        = EVENT_SM1m_delta[~np.isnan(EVENT_SM1m_delta)]
            nonevent_datax     = NONEVENT_SM1m_delta[~np.isnan(NONEVENT_SM1m_delta)]
            event_datay        = EVENT_SM1m_max[~np.isnan(EVENT_SM1m_max)]
            nonevent_datay     = NONEVENT_SM1m_max[~np.isnan(NONEVENT_SM1m_max)]
            if (len(event_datax) > 25):
                percentile_x   = calculate_bootstrap(event_datax, nonevent_datax)
                SM1m_delta_PERCENTILES[1, jdim, idim]   = percentile_x
                percentile_x   = calculate_bootstrap(event_datay, nonevent_datay)
                SM1m_anomaly_PERCENTILES[1, jdim, idim] = percentile_x
            event_datax        = EVENT_SM1mTR_delta[~np.isnan(EVENT_SM1mTR_delta)]
            nonevent_datax     = NONEVENT_SM1mTR_delta[~np.isnan(NONEVENT_SM1mTR_delta)]
            event_datay        = EVENT_SM1mTR_max[~np.isnan(EVENT_SM1mTR_max)]
            nonevent_datay     = NONEVENT_SM1mTR_max[~np.isnan(NONEVENT_SM1mTR_max)]
            if (len(event_datax) > 25):
                percentile_x   = calculate_bootstrap(event_datax, nonevent_datax)
                SM1mTR_delta_PERCENTILES[1, jdim, idim]   = percentile_x
                percentile_x   = calculate_bootstrap(event_datay, nonevent_datay)
                SM1mTR_anomaly_PERCENTILES[1, jdim, idim] = percentile_x
            event_datax        = EVENT_SM1mNTR_delta[~np.isnan(EVENT_SM1mNTR_delta)]
            nonevent_datax     = NONEVENT_SM1mNTR_delta[~np.isnan(NONEVENT_SM1mNTR_delta)]
            event_datay        = EVENT_SM1mNTR_max[~np.isnan(EVENT_SM1mNTR_max)]
            nonevent_datay     = NONEVENT_SM1mNTR_max[~np.isnan(NONEVENT_SM1mNTR_max)]
            if (len(event_datax) > 25):
                percentile_x   = calculate_bootstrap(event_datax, nonevent_datax)
                SM1mNTR_delta_PERCENTILES[1, jdim, idim]   = percentile_x
                percentile_x   = calculate_bootstrap(event_datay, nonevent_datay)
                SM1mNTR_anomaly_PERCENTILES[1, jdim, idim] = percentile_x
            
            for ievent in range(event_ind+1):
                if (EVENT_MCSmask[ievent] > 1.5):  # an MCS event
                    EVENT_SM1m_delta_copy[ievent]         = np.nan
                    EVENT_SM1m_max_copy[ievent]           = np.nan
                    NONEVENT_SM1m_delta_copy[ievent,:]    = np.nan
                    NONEVENT_SM1m_max_copy[ievent,:]      = np.nan
                    EVENT_SM1mTR_delta_copy[ievent]       = np.nan
                    EVENT_SM1mTR_max_copy[ievent]         = np.nan
                    NONEVENT_SM1mTR_delta_copy[ievent,:]  = np.nan
                    NONEVENT_SM1mTR_max_copy[ievent,:]    = np.nan
                    EVENT_SM1mNTR_delta_copy[ievent]      = np.nan
                    EVENT_SM1mNTR_max_copy[ievent]        = np.nan
                    NONEVENT_SM1mNTR_delta_copy[ievent,:] = np.nan
                    NONEVENT_SM1mNTR_max_copy[ievent,:]   = np.nan
            print('event_ind, ievent:', event_ind, ievent)
            event_datax        = EVENT_SM1m_delta_copy[~np.isnan(EVENT_SM1m_delta_copy)]
            nonevent_datax     = NONEVENT_SM1m_delta_copy[~np.isnan(NONEVENT_SM1m_delta_copy)]
            event_datay        = EVENT_SM1m_max_copy[~np.isnan(EVENT_SM1m_max_copy)]
            nonevent_datay     = NONEVENT_SM1m_max_copy[~np.isnan(NONEVENT_SM1m_max_copy)]
            if (len(event_datax) > 25):
                percentile_x   = calculate_bootstrap(event_datax, nonevent_datax)
                SM1m_delta_PERCENTILES[2, jdim, idim]   = percentile_x
                percentile_x   = calculate_bootstrap(event_datay, nonevent_datay)
                SM1m_anomaly_PERCENTILES[2, jdim, idim] = percentile_x
            event_datax        = EVENT_SM1mTR_delta_copy[~np.isnan(EVENT_SM1mTR_delta_copy)]
            nonevent_datax     = NONEVENT_SM1mTR_delta_copy[~np.isnan(NONEVENT_SM1mTR_delta_copy)]
            event_datay        = EVENT_SM1mTR_max_copy[~np.isnan(EVENT_SM1mTR_max_copy)]
            nonevent_datay     = NONEVENT_SM1mTR_max_copy[~np.isnan(NONEVENT_SM1mTR_max_copy)]
            if (len(event_datax) > 25):
                percentile_x   = calculate_bootstrap(event_datax, nonevent_datax)
                SM1mTR_delta_PERCENTILES[2, jdim, idim]   = percentile_x
                percentile_x   = calculate_bootstrap(event_datay, nonevent_datay)
                SM1mTR_anomaly_PERCENTILES[2, jdim, idim] = percentile_x
            event_datax        = EVENT_SM1mNTR_delta_copy[~np.isnan(EVENT_SM1mNTR_delta_copy)]
            nonevent_datax     = NONEVENT_SM1mNTR_delta_copy[~np.isnan(NONEVENT_SM1mNTR_delta_copy)]
            event_datay        = EVENT_SM1mNTR_max_copy[~np.isnan(EVENT_SM1mNTR_max_copy)]
            nonevent_datay     = NONEVENT_SM1mNTR_max_copy[~np.isnan(NONEVENT_SM1mNTR_max_copy)]
            if (len(event_datax) > 25):
                percentile_x   = calculate_bootstrap(event_datax, nonevent_datax)
                SM1mNTR_delta_PERCENTILES[2, jdim, idim]   = percentile_x
                percentile_x   = calculate_bootstrap(event_datay, nonevent_datay)
                SM1mNTR_anomaly_PERCENTILES[2, jdim, idim] = percentile_x


#-----save output-----
outfile0   =  path5 + 'Obtain_SMpercentiles_AfternoonRain_'+monthstr[imonth]+'_5degree1d_filterRX1.nc'
rootgrp    =  Dataset(outfile0,'w',format='NETCDF4')
p1         =  rootgrp.createDimension('time',3)
y1         =  rootgrp.createDimension('y', lat_blocks.shape[0])
x1         =  rootgrp.createDimension('x', lat_blocks.shape[1])

latitude   = rootgrp.createVariable('lat','f8',('y','x',))
latitude.long_name = 'latitude'
latitude.units     = 'degrees_north'
latitude[:,:]      = lat_blocks_center

longitude  = rootgrp.createVariable('lon','f8',('y','x',))
longitude.longname = 'longitude'
longitude.units    = 'degrees_east'
longitude[:,:]     = lon_blocks_center

nn8     = rootgrp.createVariable('SM1mdelta_percentile','f4',('time','y','x',))
nn8.long_name = 'number of event of afternoon rainfall occurs'
nn8.units     = '1'
nn8.missing_value = 1.0e+20
nn8[:,:,:]          = SM1m_delta_PERCENTILES

nn9     = rootgrp.createVariable('SM1manomaly_percentile','f4',('time','y','x',))
nn9.long_name = 'number of non-event samples of afternoon rainfall occurs'
nn9.units     = '1'
nn9.missing_value = 1.0e+20
nn9[:,:,:]          = SM1m_anomaly_PERCENTILES

pp5     = rootgrp.createVariable('SM1mTRdelta_percentile','f4',('time','y','x',))
pp5.long_name = 'percentiles of SM1 delta afternoon rainfall occurs'
pp5.units     = '1'
pp5.missing_value = 1.0e+20
pp5[:,:,:]          = SM1mTR_delta_PERCENTILES   

pp6     = rootgrp.createVariable('SM1mTRanomaly_percentile','f4',('time','y','x',))
pp6.long_name = 'percentiles of SM1 anomaly afternoon rainfall occurs'
pp6.units     = '1'
pp6.missing_value = 1.0e+20
pp6[:,:,:]          = SM1mTR_anomaly_PERCENTILES

pp7     = rootgrp.createVariable('SM1mNTRdelta_percentile','f4',('time','y','x',))
pp7.long_name = 'percentiles of SM1 delta afternoon rainfall occurs'
pp7.units     = '1'
pp7.missing_value = 1.0e+20
pp7[:,:,:]          = SM1mNTR_delta_PERCENTILES   

pp8     = rootgrp.createVariable('SM1mNTRanomaly_percentile','f4',('time','y','x',))
pp8.long_name = 'percentiles of SM1 anomaly afternoon rainfall occurs'
pp8.units     = '1'
pp8.missing_value = 1.0e+20
pp8[:,:,:]          = SM1mNTR_anomaly_PERCENTILES

rootgrp.close()
                            

    
    

