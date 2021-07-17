#!/usr/bin/env python
# coding: utf-8

# ## Script to calculate correlations between SM and ET
# 1. read in pentads SM smoothed
# 2. read in pentads ET smoothed
# 3. calculate their correlation and signficance

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

path1         = '/pic/projects/flood/CONUS_simulation/'
path5         = '/pic/projects/flood/CONUS_simulation/LDASOUT_LandAtmo/pentads/'

filename1     = path1 + 'geo_em.d01.nc'
data1         = Dataset(filename1, 'r', format='NETCDF4')
lat1          = data1.variables['CLAT']
lon1          = data1.variables['CLONG']
landmask      = data1.variables['LANDMASK']     # 0->water,  1->land
print(lon1.shape)
lon11         = lon1[0,:,:]
lat11         = lat1[0,:,:]
land11        = landmask[0,:,:]
print(lon11.shape)


# In[ ]:


#-----loop through each year to read in July data----
for imonth in range(1,5):
    ET_TOTAL_JULY     = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    ET_MCS_JULY       = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    ET_NONMCS_JULY    = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    PP_TOTAL_JULY     = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    PP_MCS_JULY       = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    PP_NONMCS_JULY    = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    SM1_TOTAL_JULY    = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    SM1_MCS_JULY      = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    SM1_NONMCS_JULY   = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    SM1m_TOTAL_JULY   = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    SM1m_MCS_JULY     = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))
    SM1m_NONMCS_JULY  = np.zeros((6*(year2-year1), lon11.shape[0], lon11.shape[1]))

    for iyear in range(year1, year2):
        filename1    = path5 + 'Obtain_ET_allpentads_'   +str(iyear)+'_160GridSmooth_corrected.nc'
        filename2    = path5 + 'Obtain_SM123_allpentads_'+str(iyear)+'_160GridSmooth.nc'
        filename3    = path5 + 'Obtain_PP_allpentads_'+str(iyear)+'_160GridSmooth.nc'
        data1        = Dataset(filename1, 'r', format='NETCDF4')
        et0          = data1.variables['ET_total_sm']
        ettr0        = data1.variables['ET_mcs_sm']  
        etntr0       = data1.variables['ET_nonmcs_sm']
        data2        = Dataset(filename2, 'r', format='NETCDF4')
        sm10         = data2.variables['SM1_total_sm']
        sm1tr0       = data2.variables['SM1_mcs_sm']
        sm1ntr0      = data2.variables['SM1_nonmcs_sm']
        sm20         = data2.variables['SM2_total_sm']
        sm2tr0       = data2.variables['SM2_mcs_sm']
        sm2ntr0      = data2.variables['SM2_nonmcs_sm']
        sm30         = data2.variables['SM3_total_sm']
        sm3tr0       = data2.variables['SM3_mcs_sm']
        sm3ntr0      = data2.variables['SM3_nonmcs_sm']
        data3        = Dataset(filename3, 'r', format='NETCDF4')
        pp0          = data3.variables['PP_total_sm']
        pptr0        = data3.variables['PP_mcs_sm']  
        ppntr0       = data3.variables['PP_nonmcs_sm']

        ind1         = (iyear-year1)*6 
        ind2         = ind1 + 6
        ind3         = imonth*6
        ind4         = ind3+6
        ET_TOTAL_JULY[ind1:ind2,:,:]    = et0[ind3:ind4,:,:]   *3600.
        ET_MCS_JULY[ind1:ind2,:,:]      = ettr0[ind3:ind4,:,:] *3600.
        ET_NONMCS_JULY[ind1:ind2,:,:]   = etntr0[ind3:ind4,:,:]*3600.
        SM1_TOTAL_JULY[ind1:ind2,:,:]   = sm10[ind3:ind4,:,:]
        SM1_MCS_JULY[ind1:ind2,:,:]     = sm1tr0[ind3:ind4,:,:]
        SM1_NONMCS_JULY[ind1:ind2,:,:]  = sm1ntr0[ind3:ind4,:,:]
        SM1m_TOTAL_JULY[ind1:ind2,:,:]  = np.add(np.add(sm10[ind3:ind4,:,:]*0.1*1000., sm20[ind3:ind4,:,:]*0.3*1000.), sm30[ind3:ind4,:,:]*0.6*1000.)
        SM1m_MCS_JULY[ind1:ind2,:,:]    = np.add(np.add(sm1tr0[ind3:ind4,:,:]*0.1*1000., sm2tr0[ind3:ind4,:,:]*0.3*1000.), sm3tr0[ind3:ind4,:,:]*0.6*1000.)
        SM1m_NONMCS_JULY[ind1:ind2,:,:] = np.add(np.add(sm1ntr0[ind3:ind4,:,:]*0.1*1000., sm2ntr0[ind3:ind4,:,:]*0.3*1000.), sm3ntr0[ind3:ind4,:,:]*0.6*1000.)
        PP_TOTAL_JULY[ind1:ind2,:,:]    = pp0[ind3:ind4,:,:]   
        PP_MCS_JULY[ind1:ind2,:,:]      = pptr0[ind3:ind4,:,:] 
        PP_NONMCS_JULY[ind1:ind2,:,:]   = ppntr0[ind3:ind4,:,:]

    print(ind2)
    
    #------calculate correlations-------
    CVALUES = np.zeros((9, lon11.shape[0], lon11.shape[1]))
    PVALUES = np.zeros((9, lon11.shape[0], lon11.shape[1]))
    CVALUES[:,:,:] = np.nan
    PVALUES[:,:,:] = np.nan
    
    CVALUES2 = np.zeros((9, lon11.shape[0], lon11.shape[1]))
    PVALUES2 = np.zeros((9, lon11.shape[0], lon11.shape[1]))
    CVALUES2[:,:,:] = np.nan
    PVALUES2[:,:,:] = np.nan
    
    CVALUES3 = np.zeros((9, lon11.shape[0], lon11.shape[1]))
    PVALUES3 = np.zeros((9, lon11.shape[0], lon11.shape[1]))
    CVALUES3[:,:,:] = np.nan
    PVALUES3[:,:,:] = np.nan
    
    I_SM1   = np.zeros((9, lon11.shape[0], lon11.shape[1]))
    I_SM1m  = np.zeros((9, lon11.shape[0], lon11.shape[1]))
    I_SM1[:,:,:]   = np.nan
    I_SM1m[:,:,:]  = np.nan
    
    I_ETPP   = np.zeros((9, lon11.shape[0], lon11.shape[1]))
    I_ETPP[:,:,:]   = np.nan

    for iy in range(0, lon11.shape[0]):
        for ix in range(0, lon11.shape[1]):
            if (land11[iy,ix] > 0.5):
                A1 = np.zeros((year2-year1)*5)     # TOTAL ET
                B1 = np.zeros((year2-year1)*5)     # MCS ET
                C1 = np.zeros((year2-year1)*5)     # NONMCS ET
                A2 = np.zeros((year2-year1)*5)     # TOTAL SURFACE SM
                B2 = np.zeros((year2-year1)*5)     # MCS SURFACE SM
                C2 = np.zeros((year2-year1)*5)     # NONMCS SURFACE SM
                A3 = np.zeros((year2-year1)*5)     # TOTAL 1M SM
                B3 = np.zeros((year2-year1)*5)     # MCS 1M SM
                C3 = np.zeros((year2-year1)*5)     # NONMCS 1M SM
                A4 = np.zeros((year2-year1)*5)     # TOTAL PP
                B4 = np.zeros((year2-year1)*5)     # MCS PP
                C4 = np.zeros((year2-year1)*5)     # NONMCS PP
                istep = -1
                
                #-----calculate historical trends------
                ET_TOTAL_SLOPE      = np.zeros(6)
                ET_MCS_SLOPE        = np.zeros(6)
                ET_NONMCS_SLOPE     = np.zeros(6)
                SM1_TOTAL_SLOPE     = np.zeros(6)
                SM1_MCS_SLOPE       = np.zeros(6)
                SM1_NONMCS_SLOPE    = np.zeros(6)
                SM1m_TOTAL_SLOPE    = np.zeros(6)
                SM1m_MCS_SLOPE      = np.zeros(6)
                SM1m_NONMCS_SLOPE   = np.zeros(6)
                PP_TOTAL_SLOPE      = np.zeros(6)
                PP_MCS_SLOPE        = np.zeros(6)
                PP_NONMCS_SLOPE     = np.zeros(6)
                for ipentad in range(6):
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),ET_TOTAL_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    ET_TOTAL_SLOPE[ipentad]                     = slope
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),ET_MCS_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    ET_MCS_SLOPE[ipentad]                       = slope
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),ET_NONMCS_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    ET_NONMCS_SLOPE[ipentad]                    = slope
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),SM1_TOTAL_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    SM1_TOTAL_SLOPE[ipentad]                    = slope
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),SM1_MCS_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    SM1_MCS_SLOPE[ipentad]                      = slope
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),SM1_NONMCS_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    SM1_NONMCS_SLOPE[ipentad]                   = slope  
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),SM1m_TOTAL_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    SM1m_TOTAL_SLOPE[ipentad]                   = slope
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),SM1m_MCS_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    SM1m_MCS_SLOPE[ipentad]                     = slope
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),SM1m_NONMCS_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    SM1m_NONMCS_SLOPE[ipentad]                  = slope
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),PP_TOTAL_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    PP_TOTAL_SLOPE[ipentad]                     = slope
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),PP_MCS_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    PP_MCS_SLOPE[ipentad]                       = slope
                    slope, intercept, r_value, p_value, std_err = stats.linregress(np.arange(year1,year2),PP_NONMCS_JULY[ipentad:6*(year2-year1):6,iy,ix])
                    PP_NONMCS_SLOPE[ipentad]                    = slope
                
                for iyear in range(year1,year2):
                    for ipentad in range(0,5):
                        istep = istep+1
                        ind1  = (iyear-year1)*6+ipentad
                        #----remove the seasonal cycle----
                        A1[istep] = ET_TOTAL_JULY[ind1, iy, ix]   - ET_TOTAL_JULY[ipentad:6*(year2-year1):6,iy,ix].mean()    - ET_TOTAL_SLOPE[ipentad]*(iyear-year1)
                        B1[istep] = ET_MCS_JULY[ind1, iy, ix]     - ET_MCS_JULY[ipentad:6*(year2-year1):6,iy,ix].mean()      - ET_MCS_SLOPE[ipentad]*(iyear-year1)
                        C1[istep] = ET_NONMCS_JULY[ind1, iy, ix]  - ET_NONMCS_JULY[ipentad:6*(year2-year1):6,iy,ix].mean()   - ET_NONMCS_SLOPE[ipentad]*(iyear-year1)
                        A2[istep] = SM1_TOTAL_JULY[ind1, iy, ix]  - SM1_TOTAL_JULY[ipentad:6*(year2-year1):6,iy,ix].mean()   - SM1_TOTAL_SLOPE[ipentad]*(iyear-year1)
                        B2[istep] = SM1_MCS_JULY[ind1, iy, ix]    - SM1_MCS_JULY[ipentad:6*(year2-year1):6,iy,ix].mean()     - SM1_MCS_SLOPE[ipentad]*(iyear-year1)
                        C2[istep] = SM1_NONMCS_JULY[ind1, iy, ix] - SM1_NONMCS_JULY[ipentad:6*(year2-year1):6,iy,ix].mean()  - SM1_NONMCS_SLOPE[ipentad]*(iyear-year1)
                        A3[istep] = SM1m_TOTAL_JULY[ind1, iy, ix] - SM1m_TOTAL_JULY[ipentad:6*(year2-year1):6,iy,ix].mean()  - SM1m_TOTAL_SLOPE[ipentad]*(iyear-year1)
                        B3[istep] = SM1m_MCS_JULY[ind1, iy, ix]   - SM1m_MCS_JULY[ipentad:6*(year2-year1):6,iy,ix].mean()    - SM1m_MCS_SLOPE[ipentad]*(iyear-year1)
                        C3[istep] = SM1m_NONMCS_JULY[ind1, iy, ix]- SM1m_NONMCS_JULY[ipentad:6*(year2-year1):6,iy,ix].mean() - SM1m_NONMCS_SLOPE[ipentad]*(iyear-year1)
                        A4[istep] = PP_TOTAL_JULY[ind1+1, iy, ix]   - PP_TOTAL_SLOPE[ipentad+1]*(iyear-year1)
                        B4[istep] = PP_MCS_JULY[ind1+1, iy, ix]     - PP_MCS_SLOPE[ipentad+1]*(iyear-year1)
                        C4[istep] = PP_NONMCS_JULY[ind1+1, iy, ix]  - PP_NONMCS_SLOPE[ipentad+1]*(iyear-year1)
                print(istep)
                #-----corelation and p values-----
                cor, pval = stats.spearmanr(A2,A1)   # PP_total ~ PP_total
                CVALUES[0,iy,ix] = cor
                PVALUES[0,iy,ix] = pval
                CVALUES[1,iy,ix] = cor
                PVALUES[1,iy,ix] = pval
                CVALUES[2,iy,ix] = cor
                PVALUES[2,iy,ix] = pval
                cor, pval = stats.spearmanr(B2,B1)   # PP_total ~ PP_MCS
                CVALUES[3,iy,ix] = cor
                PVALUES[3,iy,ix] = pval
                CVALUES[4,iy,ix] = cor
                PVALUES[4,iy,ix] = pval
                CVALUES[5,iy,ix] = cor
                PVALUES[5,iy,ix] = pval
                cor, pval = stats.spearmanr(C2,C1)   # PP_total ~ PP_nonMCS
                CVALUES[6,iy,ix] = cor
                PVALUES[6,iy,ix] = pval
                CVALUES[7,iy,ix] = cor
                PVALUES[7,iy,ix] = pval
                CVALUES[8,iy,ix] = cor
                PVALUES[8,iy,ix] = pval
                
                #-----corelation and p values for 1m SM-----
                cor, pval = stats.spearmanr(A3,A1)   # PP_total ~ PP_total
                CVALUES2[0,iy,ix] = cor
                PVALUES2[0,iy,ix] = pval
                CVALUES2[1,iy,ix] = cor
                PVALUES2[1,iy,ix] = pval
                CVALUES2[2,iy,ix] = cor
                PVALUES2[2,iy,ix] = pval
                cor, pval = stats.spearmanr(B3,B1)   # PP_total ~ PP_MCS
                CVALUES2[3,iy,ix] = cor
                PVALUES2[3,iy,ix] = pval
                CVALUES2[4,iy,ix] = cor
                PVALUES2[4,iy,ix] = pval
                CVALUES2[5,iy,ix] = cor
                PVALUES2[5,iy,ix] = pval
                cor, pval = stats.spearmanr(C3,C1)   # PP_total ~ PP_nonMCS
                CVALUES2[6,iy,ix] = cor
                PVALUES2[6,iy,ix] = pval
                CVALUES2[7,iy,ix] = cor
                PVALUES2[7,iy,ix] = pval
                CVALUES2[8,iy,ix] = cor
                PVALUES2[8,iy,ix] = pval
                
                #------correlation between ET and PP------
                cor, pval = stats.spearmanr(A4,A1)   # PP_total ~ PP_total
                CVALUES3[0,iy,ix] = cor
                PVALUES3[0,iy,ix] = pval
                cor, pval = stats.spearmanr(B4,A1)   # PP_total ~ PP_total
                CVALUES3[1,iy,ix] = cor
                PVALUES3[1,iy,ix] = pval
                cor, pval = stats.spearmanr(C4,A1)   # PP_total ~ PP_total
                CVALUES3[2,iy,ix] = cor
                PVALUES3[2,iy,ix] = pval
                cor, pval = stats.spearmanr(A4,B1)   # PP_total ~ PP_MCS
                CVALUES3[3,iy,ix] = cor
                PVALUES3[3,iy,ix] = pval
                cor, pval = stats.spearmanr(B4,B1)   # PP_total ~ PP_MCS
                CVALUES3[4,iy,ix] = cor
                PVALUES3[4,iy,ix] = pval
                cor, pval = stats.spearmanr(C4,B1)   # PP_total ~ PP_MCS
                CVALUES3[5,iy,ix] = cor
                PVALUES3[5,iy,ix] = pval
                cor, pval = stats.spearmanr(A4,C1)   # PP_total ~ PP_nonMCS
                CVALUES3[6,iy,ix] = cor
                PVALUES3[6,iy,ix] = pval
                cor, pval = stats.spearmanr(B4,C1)   # PP_total ~ PP_nonMCS
                CVALUES3[7,iy,ix] = cor
                PVALUES3[7,iy,ix] = pval
                cor, pval = stats.spearmanr(C4,C1)   # PP_total ~ PP_nonMCS
                CVALUES3[8,iy,ix] = cor
                PVALUES3[8,iy,ix] = pval
                
                
                #-----now calculate the coupling index (Dirmeyer 2011)-----
                #     I(sm) = std(sm)*(dET)/d(SM)
                slope, intercept, r_value, p_value, std_err = stats.linregress(A2,A1)
                I_SM1[0,iy,ix]   = slope * np.std(A2) 
                I_SM1[1,iy,ix]   = slope * np.std(A2)
                I_SM1[2,iy,ix]   = slope * np.std(A2)
                slope, intercept, r_value, p_value, std_err = stats.linregress(B2,B1)
                I_SM1[3,iy,ix]   = slope * np.std(B2) 
                I_SM1[4,iy,ix]   = slope * np.std(B2) 
                I_SM1[5,iy,ix]   = slope * np.std(B2) 
                slope, intercept, r_value, p_value, std_err = stats.linregress(C2,C1)
                I_SM1[6,iy,ix]   = slope * np.std(C2) 
                I_SM1[7,iy,ix]   = slope * np.std(C2) 
                I_SM1[8,iy,ix]   = slope * np.std(C2) 
                
                slope, intercept, r_value, p_value, std_err = stats.linregress(A3,A1)
                I_SM1m[0,iy,ix]   = slope * np.std(A3) 
                I_SM1m[1,iy,ix]   = slope * np.std(A3) 
                I_SM1m[2,iy,ix]   = slope * np.std(A3)
                slope, intercept, r_value, p_value, std_err = stats.linregress(B3,B1)
                I_SM1m[3,iy,ix]   = slope * np.std(B3) 
                I_SM1m[4,iy,ix]   = slope * np.std(B3) 
                I_SM1m[5,iy,ix]   = slope * np.std(B3) 
                slope, intercept, r_value, p_value, std_err = stats.linregress(C3,C1)
                I_SM1m[6,iy,ix]   = slope * np.std(C3) 
                I_SM1m[7,iy,ix]   = slope * np.std(C3)
                I_SM1m[8,iy,ix]   = slope * np.std(C3) 
                
                slope, intercept, r_value, p_value, std_err = stats.linregress(A1,A4)
                I_ETPP[0,iy,ix]   = slope  
                slope, intercept, r_value, p_value, std_err = stats.linregress(A1,B4)
                I_ETPP[1,iy,ix]   = slope  
                slope, intercept, r_value, p_value, std_err = stats.linregress(A1,C4)
                I_ETPP[2,iy,ix]   = slope
                slope, intercept, r_value, p_value, std_err = stats.linregress(B1,A4)
                I_ETPP[3,iy,ix]   = slope
                slope, intercept, r_value, p_value, std_err = stats.linregress(B1,B4)
                I_ETPP[4,iy,ix]   = slope 
                slope, intercept, r_value, p_value, std_err = stats.linregress(B1,C4)
                I_ETPP[5,iy,ix]   = slope
                slope, intercept, r_value, p_value, std_err = stats.linregress(C1,A4)
                I_ETPP[6,iy,ix]   = slope  
                slope, intercept, r_value, p_value, std_err = stats.linregress(C1,B4)
                I_ETPP[7,iy,ix]   = slope 
                slope, intercept, r_value, p_value, std_err = stats.linregress(C1,C4)
                I_ETPP[8,iy,ix]   = slope
                
    #-----multiply SM->ET, ET->PP-----
    I_SM1PP  = np.multiply(I_SM1,  I_ETPP)
    I_SM1mPP = np.multiply(I_SM1m, I_ETPP)
    CVALUES_SM1PP  = np.multiply(CVALUES,  CVALUES3)
    CVALUES_SM1mPP = np.multiply(CVALUES2, CVALUES3)
                
    #------output these matrices as netcdf files--------

    outfile0   =  path1 + '/LDASOUT_LandAtmo/Obtain_'+ monthstr[imonth] +'SM5ET5XET5PP5_correlations_160GridSmooth_detrend.nc'
    rootgrp    =  Dataset(outfile0,'w',format='NETCDF4')
    p1         =  rootgrp.createDimension('pentad',9)
    y1         =  rootgrp.createDimension('y', lon11.shape[0])
    x1         =  rootgrp.createDimension('x', lon11.shape[1])

    latitude   = rootgrp.createVariable('lat','f8',('y','x',))
    latitude.long_name = 'latitude'
    latitude.units     = 'degrees_north'
    latitude[:,:]      = lat11

    longitude  = rootgrp.createVariable('lon','f8',('y','x',))
    longitude.longname = 'longitude'
    longitude.units    = 'degrees_east'
    longitude[:,:]     = lon11

    cor1   = rootgrp.createVariable('CVALUE_SM1ET','f4',('pentad','y','x',))
    cor1.long_name = 'Correlations for lagged precip'
    cor1.units     = '1'
    cor1.missing_value = 1.0e+20
    cor1[:,:,:]        = CVALUES

    pvar1   = rootgrp.createVariable('PVALUE_SM1ET','f4',('pentad','y','x',))
    pvar1.long_name = 'P values for lagged precip'
    pvar1.units     = '1'
    pvar1.missing_value = 1.0e+20
    pvar1[:,:,:]        = PVALUES
    
    cor2   = rootgrp.createVariable('CVALUE_SM1mET','f4',('pentad','y','x',))
    cor2.long_name = 'Correlations for lagged precip'
    cor2.units     = '1'
    cor2.missing_value = 1.0e+20
    cor2[:,:,:]        = CVALUES2

    pvar2   = rootgrp.createVariable('PVALUE_SM1mET','f4',('pentad','y','x',))
    pvar2.long_name = 'P values for lagged precip'
    pvar2.units     = '1'
    pvar2.missing_value = 1.0e+20
    pvar2[:,:,:]        = PVALUES2
    
    cor3   = rootgrp.createVariable('CVALUE_ETPP','f4',('pentad','y','x',))
    cor3.long_name = 'Correlations for lagged precip'
    cor3.units     = '1'
    cor3.missing_value = 1.0e+20
    cor3[:,:,:]        = CVALUES3
    
    cor4   = rootgrp.createVariable('CVALUE_SM1PP','f4',('pentad','y','x',))
    cor4.long_name = 'Correlations for lagged precip'
    cor4.units     = '1'
    cor4.missing_value = 1.0e+20
    cor4[:,:,:]        = CVALUES_SM1PP
    
    cor5   = rootgrp.createVariable('CVALUE_SM1mPP','f4',('pentad','y','x',))
    cor5.long_name = 'Correlations for lagged precip'
    cor5.units     = '1'
    cor5.missing_value = 1.0e+20
    cor5[:,:,:]        = CVALUES_SM1mPP
    
    Ivar1   = rootgrp.createVariable('I_SM1ET','f4',('pentad','y','x',))
    Ivar1.long_name = 'coupling index for surface SM and ET'
    Ivar1.units     = '1'
    Ivar1.missing_value = 1.0e+20
    Ivar1[:,:,:]        = I_SM1
    
    Ivar2   = rootgrp.createVariable('I_SM1mET','f4',('pentad','y','x',))
    Ivar2.long_name = 'coupling index for 1m SM and ET'
    Ivar2.units     = '1'
    Ivar2.missing_value = 1.0e+20
    Ivar2[:,:,:]        = I_SM1m
    
    Ivar3   = rootgrp.createVariable('I_ETPP','f4',('pentad','y','x',))
    Ivar3.long_name = 'coupling index for 1m SM and ET'
    Ivar3.units     = '1'
    Ivar3.missing_value = 1.0e+20
    Ivar3[:,:,:]        = I_ETPP
    
    Ivar4   = rootgrp.createVariable('I_SM1PP','f4',('pentad','y','x',))
    Ivar4.long_name = 'coupling index for 1m SM and ET'
    Ivar4.units     = '1'
    Ivar4.missing_value = 1.0e+20
    Ivar4[:,:,:]        = I_SM1PP
    
    Ivar5   = rootgrp.createVariable('I_SM1mPP','f4',('pentad','y','x',))
    Ivar5.long_name = 'coupling index for 1m SM and ET'
    Ivar5.units     = '1'
    Ivar5.missing_value = 1.0e+20
    Ivar5[:,:,:]        = I_SM1mPP

    rootgrp.close()
    print('Output saved.')

