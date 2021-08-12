#!/usr/bin/env python
# coding: utf-8

# ## Script to obtain afternoon rainfall occurence and its soil moisture conditions
# 1. divide the CONUS domain into 120km*120km blocks (30*30 pixels)
# 2. loop through each box and each day
# 3. accumulate afternoon rainfall and smooth into 24km*24km boxes
# 4. find the box with maximum afternoon rainfall if > 2mm
# 5. adjust the block according to the location of maximum box
# 6. find the box/multiple-box with minimum afternoon rainfall
# 7. record their location, SM
# 

# In[ ]:


import os
import sys
import numpy as np
from netCDF4 import Dataset
import csv
import string
from astropy.convolution import Box2DKernel, convolve

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
height11  = height[0,:,:]
print(lon11.shape)


# In[ ]:


#----do smooth for lat and lon matrices----
box_size       = 6   #  4km*6 =  24km
block_size     = 5   # 24km*5 = 120km

lon_boxes      = np.zeros((int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))
lat_boxes      = np.zeros((int(lon11.shape[0]/box_size), int(lon11.shape[1]/box_size)))

lon_1          = np.mean(np.reshape(lon11, (    lon11.shape[0],                     int(lon11.shape[1]/box_size), box_size)), axis=2)
lon_2          = np.mean(np.reshape(lon_1, (int(lon11.shape[0]/box_size), box_size, int(lon11.shape[1]/box_size)          )), axis=1)
print(lon_2.shape)

lat_1          = np.mean(np.reshape(lat11, (    lon11.shape[0],                     int(lon11.shape[1]/box_size), box_size)), axis=2)
lat_2          = np.mean(np.reshape(lat_1, (int(lon11.shape[0]/box_size), box_size, int(lon11.shape[1]/box_size)          )), axis=1)

#----also process height and water body-----
height_1       = np.mean(np.reshape(height11,(    lon11.shape[0],                     int(lon11.shape[1]/box_size), box_size)), axis=2)
height_2       = np.mean(np.reshape(height_1,(int(lon11.shape[0]/box_size), box_size, int(lon11.shape[1]/box_size)          )), axis=1)

land_1         = np.mean(np.reshape(land11, (    lon11.shape[0],                     int(lon11.shape[1]/box_size), box_size)), axis=2)
land_2         = np.mean(np.reshape(land_1, (int(lon11.shape[0]/box_size), box_size, int(lon11.shape[1]/box_size)          )), axis=1)

jind_matrix    = np.zeros(lat_2.shape)
iind_matrix    = np.zeros(lat_2.shape)
for jj in range(jind_matrix.shape[0]):
    jind_matrix[jj,:] = jj
for ii in range(iind_matrix.shape[1]):
    iind_matrix[:,ii] = ii


# In[ ]:


import datetime

def find_timestring(y1,m1,d1,h1):
    time1 = datetime.datetime(y1,m1,d1)
    time2 = time1 + datetime.timedelta(hours=h1)
    return time2.strftime('%Y%m%d%H')


# In[ ]:


afternoon_rain                = 3   # set the criteria > 2 mm

frequency_matrix_month        = np.zeros((5, lon_2.shape[0], lon_2.shape[1]))
frequency_mcs_month           = np.zeros((5, lon_2.shape[0], lon_2.shape[1]))

for iyear in range(year1, year2):
    MAX_IND_MATRIX            = np.zeros((sum(days), lon_2.shape[0], lon_2.shape[1]))
    MAX_SUR_MATRIX            = np.zeros((sum(days), lon_2.shape[0], lon_2.shape[1]))
    DAILYPP_MATRIX            = np.zeros((sum(days), lon_2.shape[0], lon_2.shape[1]))
    MAX_MCSMASK_MATRIX        = np.zeros((sum(days), lon_2.shape[0], lon_2.shape[1]))
    MORNINGMASK_MATRIX        = np.zeros((sum(days), lon_2.shape[0], lon_2.shape[1]))
    MAX_IND_MATRIX[:,:,:]     = np.nan 
    MAX_SUR_MATRIX[:,:,:]     = np.nan
    day_ind                   = -1
    
    path2     = path1 + 'LDASOUT_' + str(iyear)
    for imonth in range(5):
        max_ind              = -1
        TIME_LIST            = []
        MAX_LATLOC_LIST      = []
        MAX_LONLOC_LIST      = []
        MIN_LATLOC_LIST      = []
        MIN_LONLOC_LIST      = []
        remove_list          = []
        pp_max_list      = []
        for iday in range(1, days[imonth]+1):

            day_ind          = day_ind+1
            hour_start       = 27       # local time: 27-6=21
            hour_end         = 33       # local time: 33-6=27
            rain_6hours      = np.zeros((lon11.shape))
            rain_6hours_mcs  = np.zeros((lon11.shape))
            rain_morning     = np.zeros((lon11.shape))
            for ihour in range(hour_start, hour_end):
                datetime_str = find_timestring(iyear,imonth+4,iday,ihour)
                filename1    = path2 + '/' + datetime_str + '.LDASOUT_DOMAIN1.nc'
                print(filename1)
                data1        = Dataset(filename1,'r',format='NETCDF4')
                precip       = data1.variables['RAINRATE']
                precip_tr    = data1.variables['RAINRATE_TR']
                rain_6hours  = np.add(rain_6hours, precip[0,:,:])
                rain_6hours_mcs = np.add(rain_6hours_mcs, precip_tr[0,:,:])
            #----also calculate morning rainfall-----
            hour_start1      = 15      # local time: 15-6= 9
            hour_end1        = 27      # local time: 27-6= 21   #afternoon hours
            datetime_str_am  = find_timestring(iyear,imonth+4,iday,hour_start1)
            for ihour in range(hour_start1, hour_end1):
                datetime_str = find_timestring(iyear,imonth+4,iday,ihour)
                filename1    = path2 + '/' + datetime_str + '.LDASOUT_DOMAIN1.nc'
                print('Morning rain files:',filename1)
                data1        = Dataset(filename1,'r',format='NETCDF4')
                precip       = data1.variables['RAINRATE']
                rain_morning = np.add(rain_morning, precip[0,:,:])
            #----smooth into 24-km by 24-km grids-----
            rain_6hours      = np.where(land11>0.5, rain_6hours, np.nan)
            rain_1           = np.nanmean(np.reshape(rain_6hours,(lon11.shape[0],              int(lon11.shape[1]/box_size),                  box_size)),axis=2)
            rain_2           = np.nanmean(np.reshape(rain_1,     (int(lon11.shape[0]/box_size),                   box_size,    lon11.shape[1]/box_size)),axis=1)
            print(rain_2.shape)
            rain_6hours_mcs  = np.where(land11>0.5, rain_6hours_mcs, np.nan)
            rain_1_mcs       = np.nanmean(np.reshape(rain_6hours_mcs,(lon11.shape[0],              int(lon11.shape[1]/box_size),                  box_size)),axis=2)
            rain_2_mcs       = np.nanmean(np.reshape(rain_1_mcs,     (int(lon11.shape[0]/box_size),                   box_size,    lon11.shape[1]/box_size)),axis=1)
            rain_morning     = np.where(land11>0.5, rain_morning, np.nan)
            rain_1_morning   = np.nanmean(np.reshape(rain_morning,       (lon11.shape[0],              int(lon11.shape[1]/box_size),                  box_size)),axis=2)
            rain_2_morning   = np.nanmean(np.reshape(rain_1_morning,     (int(lon11.shape[0]/box_size),                   box_size,    lon11.shape[1]/box_size)),axis=1)
            
            DAILYPP_MATRIX[day_ind,:,:]     = rain_2
            MORNINGMASK_MATRIX[day_ind,:,:] = np.where(rain_2_morning > 2.0, 2, 1)
            
            for jbox in range(0, int(rain_2.shape[0])):
                for ibox in range(0, int(rain_2.shape[1])):
                    pp_box           = rain_2[jbox, ibox]
                    pp_block         = rain_2[max(0,jbox-int(block_size/2)):min(rain_2.shape[0],jbox+int(block_size/2)+1),max(0,ibox-int(block_size/2)):min(rain_2.shape[1],ibox+int(block_size/2)+1)]
                    j_block          = jind_matrix[max(0,jbox-int(block_size/2)):min(rain_2.shape[0],jbox+int(block_size/2)+1),max(0,ibox-int(block_size/2)):min(rain_2.shape[1],ibox+int(block_size/2)+1)]
                    i_block          = iind_matrix[max(0,jbox-int(block_size/2)):min(rain_2.shape[0],jbox+int(block_size/2)+1),max(0,ibox-int(block_size/2)):min(rain_2.shape[1],ibox+int(block_size/2)+1)]
                    if (pp_box > afternoon_rain):
                        jloc_old     = int(jbox)
                        iloc_old     = int(ibox)
                        j_block_old  = np.copy(j_block)
                        i_block_old  = np.copy(i_block)
                        pp_block_old = np.copy(pp_block)
                        #-----find the maximum in the block and adjust block location-----
                        max_loc      = np.unravel_index(np.nanargmax(pp_block_old),pp_block_old.shape)
                        jloc_new     = int(j_block_old[max_loc[0], max_loc[1]])
                        iloc_new     = int(i_block_old[max_loc[0], max_loc[1]])
                        pp_block_new = np.copy(pp_block_old)
                        j_block_new  = np.copy(j_block_old)
                        i_block_new  = np.copy(i_block_old)
                        while (jloc_new != jloc_old or iloc_new != iloc_old):    # max is not at the center, iterate
                            jloc_old     = jloc_new
                            iloc_old     = iloc_new
                            pp_block_new = rain_2[max(0,jloc_new-int(block_size/2)):min(rain_2.shape[0],jloc_new+int(block_size/2)+1),max(0,iloc_new-int(block_size/2)):min(rain_2.shape[1],iloc_new+int(block_size/2)+1)]
                            j_block_new  = jind_matrix[max(0,jloc_new-int(block_size/2)):min(rain_2.shape[0],jloc_new+int(block_size/2)+1),max(0,iloc_new-int(block_size/2)):min(rain_2.shape[1],iloc_new+int(block_size/2)+1)]
                            i_block_new  = iind_matrix[max(0,jloc_new-int(block_size/2)):min(rain_2.shape[0],jloc_new+int(block_size/2)+1),max(0,iloc_new-int(block_size/2)):min(rain_2.shape[1],iloc_new+int(block_size/2)+1)]
                            max_loc      = np.unravel_index(np.nanargmax(pp_block_new),pp_block_new.shape)
                            jloc_new     = int(j_block_new[max_loc[0], max_loc[1]])
                            iloc_new     = int(i_block_new[max_loc[0], max_loc[1]])
                            
                        print('MAX_IND_MAXTRIX: ', MAX_IND_MATRIX[day_ind, jloc_new, iloc_new], 'day_ind', day_ind)
                        #-----tell if the rain is associated with MCS or not----
                        ppmcs_box        = rain_2_mcs[jloc_new, iloc_new]
                        mcs_mask         = 1      # no mcs
#                         if (ppmcs_box > 1.0):
                        if (ppmcs_box > rain_2[jloc_new, iloc_new]*0.5):
                            mcs_mask     = 2       # yes mcs
                        else:
                            mcs_mask     = 1
                        #-----tell if there is morning rainfall in this block-----
                        pp_morning_block = rain_2_morning[max(0,jloc_new-int(block_size/2)):min(rain_2.shape[0],jloc_new+int(block_size/2)+1),max(0,iloc_new-int(block_size/2)):min(rain_2.shape[1],iloc_new+int(block_size/2)+1)]
                        morning_mask     = 1      # no moring rain
                        if (pp_morning_block.max() > 2.0):
                            morning_mask = 2      # yes morning rain
                        else:
                            morning_mask = 1
                        #-----tell if there is topography-----
                        height_block     = height_2[max(0,jloc_new-int(block_size/2)):min(rain_2.shape[0],jloc_new+int(block_size/2)+1),max(0,iloc_new-int(block_size/2)):min(rain_2.shape[1],iloc_new+int(block_size/2)+1)]
                        dims_1           = height_block.shape[0]
                        dims_2           = height_block.shape[1]
                        topo_mask        = 1
                        if ((height_block.max()-height_block.min()) > 300./block_size*min(dims_1, dims_2)):
                            topo_mask    = 2     # yes with terrain
                        else:
                            topo_mask    = 1
                        #-----tell if there is water body in the block-----
                        land_block       = land_2[max(0,jloc_new-int(block_size/2)):min(rain_2.shape[0],jloc_new+int(block_size/2)+1),max(0,iloc_new-int(block_size/2)):min(rain_2.shape[1],iloc_new+int(block_size/2)+1)]
                        land_mask        = 1
                        if (land_block.min() < 0.95):
                            land_mask    = 2
                        else:
                            land_mask    = 1
                            
                        #-----filter out all these: morning rainfall, topo, water body----
                        if ( topo_mask < 1.5 and land_mask<1.5):
                            #----tell if the block is all nan-----
                            event_block = MAX_SUR_MATRIX[day_ind, max(0,jloc_new-int(block_size/2)):min(rain_2.shape[0],jloc_new+int(block_size/2)+1),max(0,iloc_new-int(block_size/2)):min(rain_2.shape[1],iloc_new+int(block_size/2)+1)]
                            if (len(event_block[~np.isnan(event_block)])>0.5):
                                #----some of the block has been occupied, meaning there is overlap with another event----
                                print('overlap event detected.')
                                #----only keep the larger event----
                                for event_ind in np.unique(event_block[~np.isnan(event_block)]):
                                    print('event_ind, max_ind, len(pp_max_list):', event_ind)
                                    pp_max_previous = pp_max_list[int(event_ind)]
                                    if (pp_max_previous < np.nanargmax(pp_block_new)):
                                        #----replace this block-----
                                        print('replacing previous event: pp_old, pp_new:', pp_max_previous, np.nanargmax(pp_block_new))
                                        MAX_IND_MATRIX[day_ind,:,:] = np.where(MAX_IND_MATRIX[day_ind,:,:] == event_ind, np.nan, MAX_IND_MATRIX[day_ind,:,:])
                                        MAX_SUR_MATRIX[day_ind,:,:] = np.where(MAX_SUR_MATRIX[day_ind,:,:] == event_ind, np.nan, MAX_SUR_MATRIX[day_ind,:,:])
                                        max_ind      = max_ind+1
                                        MAX_IND_MATRIX[day_ind, jloc_new, iloc_new] = max_ind
                                        MAX_SUR_MATRIX[day_ind,max(0,jloc_new-int(block_size/2)):min(rain_2.shape[0],jloc_new+int(block_size/2)+1),max(0,iloc_new-int(block_size/2)):min(rain_2.shape[1],iloc_new+int(block_size/2)+1)] = max_ind
                                        pp_max_list.append(np.nanargmax(pp_block_new))
                                        remove_list.append(int(event_ind))
                                        MAX_LATLOC_LIST.append(jloc_new)
                                        MAX_LONLOC_LIST.append(iloc_new)
                                        TIME_LIST.append(datetime_str_am)
                                        MAX_MCSMASK_MATRIX[day_ind, jloc_new, iloc_new] = mcs_mask
                                        print('MAX_IND_MAXTRIX: ', MAX_IND_MATRIX[day_ind, jloc_new, iloc_new])
                                        pp_min       = np.nanmin(pp_block_new)
                                        min_loc      = np.where(pp_block_new == pp_min)
                                        if len(min_loc[0])>1:
                                            i_list = []
                                            j_list = []
                                            for imin in range(len(min_loc[0])):
                                                i_list.append(i_block_new[min_loc[0][imin],min_loc[1][imin]])
                                                j_list.append(j_block_new[min_loc[0][imin],min_loc[1][imin]])
                                        else:
                                            i_list = i_block_new[min_loc[0],min_loc[1]]
                                            j_list = j_block_new[min_loc[0],min_loc[1]]
                                        print(j_list)
                                        print(i_list)
                                        MIN_LATLOC_LIST.append(j_list)
                                        MIN_LONLOC_LIST.append(i_list)
                                    else:
                                        #----does not need to replace the previous event-----
                                        print('does not need to replace previous event:pp_old, pp_new:', pp_max_previous, np.nanargmax(pp_block_new))
                            else:     #all fresh grids
                                max_ind      = max_ind+1
                                MAX_IND_MATRIX[day_ind, jloc_new, iloc_new] = max_ind
                                MAX_SUR_MATRIX[day_ind,max(0,jloc_new-int(block_size/2)):min(rain_2.shape[0],jloc_new+int(block_size/2)+1),max(0,iloc_new-int(block_size/2)):min(rain_2.shape[1],iloc_new+int(block_size/2)+1)] = max_ind
                                pp_max_list.append(np.nanargmax(pp_block_new))
                                MAX_LATLOC_LIST.append(jloc_new)
                                MAX_LONLOC_LIST.append(iloc_new)
                                TIME_LIST.append(datetime_str_am)
                                MAX_MCSMASK_MATRIX[day_ind, jloc_new, iloc_new] = mcs_mask
                                MORNINGMASK_MATRIX[day_ind, jloc_new, iloc_new] = morning_mask
                                print('new record box:')
                                print(int(MAX_LATLOC_LIST[-1]), int(MAX_LONLOC_LIST[-1]))
                                pp_min       = np.nanmin(pp_block_new)
                                min_loc      = np.where(pp_block_new == pp_min)

                                print('MAX_IND_MAXTRIX: ', MAX_IND_MATRIX[day_ind, jloc_new, iloc_new])
                                if len(min_loc[0])>1:
                                    i_list = []
                                    j_list = []
                                    for imin in range(len(min_loc[0])):
                                        i_list.append(i_block_new[min_loc[0][imin],min_loc[1][imin]])
                                        j_list.append(j_block_new[min_loc[0][imin],min_loc[1][imin]])
                                else:
                                    i_list = i_block_new[min_loc[0],min_loc[1]]
                                    j_list = j_block_new[min_loc[0],min_loc[1]]
                                print(j_list)
                                print(i_list)
                                MIN_LATLOC_LIST.append(j_list)
                                MIN_LONLOC_LIST.append(i_list)
                        else:
                            print('Block excluded. morning_mask, topo_mask, land_mask:',morning_mask, topo_mask, land_mask)

        #----remove replaced blocks due to overlapping-----
        print('remove_list:', remove_list)
        for remove_ind in remove_list:
            del MAX_LATLOC_LIST[remove_ind]
            del MAX_LONLOC_LIST[remove_ind]
            del MIN_LATLOC_LIST[remove_ind]
            del MIN_LONLOC_LIST[remove_ind]
            del TIME_LIST[remove_ind]
            
        #----output lists for each month-----
        outfile1 = path5 + 'Night_PP_max_latind_'+ str(iyear) + months[imonth] + '_24km120km_iterate_filterRX1.txt'
        outfile2 = path5 + 'Night_PP_max_lonind_'+ str(iyear) + months[imonth] + '_24km120km_iterate_filterRX1.txt'
        outfile3 = path5 + 'Night_PP_min_latind_'+ str(iyear) + months[imonth] + '_24km120km_iterate_filterRX1.txt'
        outfile4 = path5 + 'Night_PP_min_lonind_'+ str(iyear) + months[imonth] + '_24km120km_iterate_filterRX1.txt'
        outfile5 = path5 + 'Night_PP_datetime_'  + str(iyear) + months[imonth] + '_24km120km_iterate_filterRX1.txt'
        f = open(outfile1, 'w')
        f.write('\n'.join('%s' % x for x in MAX_LATLOC_LIST))
        f.close()
        f = open(outfile2, 'w')
        f.write('\n'.join('%s' % x for x in MAX_LONLOC_LIST))
        f.close()
        f = open(outfile3, 'w')
        for k in range(len(MIN_LATLOC_LIST)):
            f.write(','.join('%s' % x for x in MIN_LATLOC_LIST[k]))
            f.write('\n')
        f.close()
        f = open(outfile4, 'w')
        for k in range(len(MIN_LATLOC_LIST)):
            f.write(','.join('%s' % x for x in MIN_LONLOC_LIST[k]))
            f.write('\n')
        f.close()
        f = open(outfile5, 'w')
        f.write('\n'.join('%s' % x for x in TIME_LIST))
        f.close()
    
    #-----output matrix as netcdf file-----
    outfile0   =  path5 + 'Obtain_NightPP_maxind_'+str(iyear)+'_24km120km_iterate_filterRX1.nc'
    rootgrp    =  Dataset(outfile0,'w',format='NETCDF4')
    p1         =  rootgrp.createDimension('time',sum(days))
    y1         =  rootgrp.createDimension('y', lon_2.shape[0])
    x1         =  rootgrp.createDimension('x', lon_2.shape[1])
    
    latitude   = rootgrp.createVariable('lat','f8',('y','x',))
    latitude.long_name = 'latitude'
    latitude.units     = 'degrees_north'
    latitude[:,:]      = lat_2

    longitude  = rootgrp.createVariable('lon','f8',('y','x',))
    longitude.longname = 'longitude'
    longitude.units    = 'degrees_east'
    longitude[:,:]     = lon_2

    pp1     = rootgrp.createVariable('MAXPP_ind','f4',('time','y','x',))
    pp1.long_name = 'index indicating where maximum afternoon rainfall occurs'
    pp1.units     = '1'
    pp1.missing_value = 1.0e+20
    pp1[:,:,:]        = MAX_IND_MATRIX
    
    pp2     = rootgrp.createVariable('PP_amount','f4',('time','y','x',))
    pp2.long_name = 'amount of afternoon rainfall'
    pp2.units     = 'mm'
    pp2.missing_value = 1.0e+20
    pp2[:,:,:]        = DAILYPP_MATRIX
    
    pp3     = rootgrp.createVariable('MAX_mcsmask','f4',('time','y','x',))
    pp3.long_name = 'index indicating if maximum rainfall is associate with MCS (1) or not (0)'
    pp3.units     = '1'
    pp3.missing_value = 1.0e+20
    pp3[:,:,:]        = MAX_MCSMASK_MATRIX
    
    pp4     = rootgrp.createVariable('MAX_morningmask','f4',('time','y','x',))
    pp4.long_name = 'index indicating if maximum rainfall is associate with morning rainfall (1) or not (0)'
    pp4.units     = '1'
    pp4.missing_value = 1.0e+20
    pp4[:,:,:]        = MORNINGMASK_MATRIX
    
    pp5     = rootgrp.createVariable('MAXPP_sur','f4',('time','y','x',))
    pp5.long_name = 'index indicating which rain event each rainy grid is associated with'
    pp5.units     = '1'
    pp5.missing_value = 1.0e+20
    pp5[:,:,:]        = MAX_SUR_MATRIX
    
    rootgrp.close()
    print('Output saved for year '+str(iyear))
    
    #-----accumulate the frequency of rain events and associated mcs in each month
    for imonth in range(5):
        ind1 = sum(days[0:imonth])
        ind2 = sum(days[0:imonth+1])
        print('imonth, ind1, ind2:', imonth, ind1,ind2)
        ones_matrix     = np.where(MAX_IND_MATRIX[ind1:ind2,:,:]>0.5, 1, 0)
        rain_days       = np.sum(ones_matrix, axis=0)
        ones_matrix_mcs = np.where(MAX_MCSMASK_MATRIX[ind1:ind2,:,:]>1.5, 1, 0)
        mcs_days        = np.sum(ones_matrix_mcs, axis=0)
        frequency_matrix_month[imonth,:,:]  = np.add(frequency_matrix_month[imonth,:,:], rain_days)
        frequency_mcs_month[imonth,:,:]     = np.add(frequency_mcs_month[imonth,:,:],  mcs_days)
        

frequency_matrix_month = np.where(frequency_matrix_month==0, np.nan, frequency_matrix_month)
frequency_mcs_month    = np.where(frequency_matrix_month==0, np.nan, frequency_mcs_month)

#-----output the rainy frequency and non-rainy frequncy matrices-----
outfile0   =  path5 + 'Obtain_NightPP_frequency_allyears_24km120km_iterate_filterRX1.nc'
rootgrp    =  Dataset(outfile0,'w',format='NETCDF4')
p1         =  rootgrp.createDimension('time',5)
y1         =  rootgrp.createDimension('y', lon_2.shape[0])
x1         =  rootgrp.createDimension('x', lon_2.shape[1])

latitude   = rootgrp.createVariable('lat','f8',('y','x',))
latitude.long_name = 'latitude'
latitude.units     = 'degrees_north'
latitude[:,:]      = lat_2

longitude  = rootgrp.createVariable('lon','f8',('y','x',))
longitude.longname = 'longitude'
longitude.units    = 'degrees_east'
longitude[:,:]     = lon_2

pp1     = rootgrp.createVariable('frequency_rain','f4',('time','y','x',))
pp1.long_name = 'rainy frequency'
pp1.units     = 'days'
pp1.missing_value = 1.0e+20
pp1[:,:,:]        = frequency_matrix_month

pp2     = rootgrp.createVariable('frequency_mcs','f4',('time','y','x',))
pp2.long_name = 'rainy frequency associated with mcs or not'
pp2.units     = 'days'
pp2.missing_value = 1.0e+20
pp2[:,:,:]        = frequency_mcs_month

rootgrp.close()

