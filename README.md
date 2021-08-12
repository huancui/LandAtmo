## The role of MCS and non-MCS rainfall in soil moisture-precipitation feedbacks
This repository includes key scripts to analyze sub-daily and pentad interactions between soil moisture and rainfall from Mesoscale Convective Systems(MCSs) and non-MCS components.

Reference: To be added.

------
## List of scripts:

- [1. ObtainSH_Locations_for_AfternoonRain_24km120km_iterate.py](#head1)
- [2. Obtain_StormEvent_Episode_MCS.py](#head2)
- [3. Obtain_StormEvent_Episode_nonMCS.py](#head3)
- [Contact us](#head4)

---------
### <a name="head1"></a>1. ObtainSH_Locations_for_AfternoonRain_24km120km_iterate.py

The script regrids the original 4-km grids into 24-km grids and looks for locations with afternoon rainfall (3pm-9pm local time). The locations with maximum and minimum rainfall within each 120-km by 120-km box are recorded in order to analyze the spatial gradients in the next step. Note that mountainous and water grids are filtered out.    

#### Input files:
|  | File names | Description |
| ----- | ------ | ------ |
| A	| %Y%m%d%H.LDASOUT_DOMAIN1.nc | read in precipitation that has been seperated into MCS and non-MCS components |
| B | geo_em.d01.nc | read in the height and land/water mask in order to filter out mountainous and water grids|

#### Output files:
|  | File names | Description |
| ----- | ------ | ------ |
| A	| Afternoon_PP_max_latind_yyyymm_24km120km_iterate_filterRX1.txt | record the latitude indices for locations with maximum rainfall in each 120kmx120km box|
| B | Afternoon_PP_max_lonind_yyyymm_24km120km_iterate_filterRX1.txt | record the longitude indices for locations with maximum rainfall in each 120kmx120km box|
| C	| Afternoon_PP_min_latind_yyyymm_24km120km_iterate_filterRX1.txt | record the latitude indices for locations with minimum rainfall in each 120kmx120km box|
| D | Afternoon_PP_min_lonind_yyyymm_24km120km_iterate_filterRX1.txt | record the longitude indices for locations with minimum rainfall in each 120kmx120km box|
| E | Afternoon_datetime_yyyymm_24km120km_iterate_filterRX1.txt | record the date and time when maximum rainfall occurs |
| F	| Obtain_AfternoonPP_maxind_yyyy_24km120km_iterate_filterRX1.nc | map the locations with maximum rainfall|
| G | Obtain_AfternoonPP_frequency_allyears_24km120km_iterate_filterRX1.nc | aggregate spatially for frequency of maximum rainfall in all years |

---------
### <a name="head2"></a>2. ObtainSH_Locations_for_NightRain_24km120km_iterate.py

The script is similar to "ObtainSH_Locations_for_AfternoonRain_24km120km_iterate.py" but records locations with nighttime rainfall (9pm-3am). 

---------
### <a name="head3"></a>3. ObtainSH_SMpercentile_AfternoonRain_5degree1d.py

The script calculates the spatial and temporal preferences when afternoon rainfall occurs. The results are saved as percentiles of temporal/spatial soil moisture anomaly when afternoon rainfall occurs in each 5-degree box. To better reflect the spatial variation, results for 5-degree boxes centered at each 1-degree grid are calculated.

#### Input files:
|  | File names | Description |
| ----- | ------ | ------ |
| A	| Obtain_9amSM1m_yyyy_24km120km.nc | read in morning soil moisture that has the MCS and non-MCS components |
| B | Obtain_AfternoonPP_maxind_yyyy_24km120km_iterate_filterRX1.nc | read in locations and rainfall amount of maximum rainfall |
| C	| Afternoon_PP_max_latind_yyyymm_24km120km_iterate_filterRX1.txt | read in latitude indices for locations with maximum rainfall in each 120kmx120km box|
| D | Afternoon_PP_max_lonind_yyyymm_24km120km_iterate_filterRX1.txt | read in longitude indices for locations with maximum rainfall in each 120kmx120km box|
| E	| Afternoon_PP_min_latind_yyyymm_24km120km_iterate_filterRX1.txt | read in latitude indices for locations with minimum rainfall in each 120kmx120km box|
| F | Afternoon_PP_min_lonind_yyyymm_24km120km_iterate_filterRX1.txt | read in longitude indices for locations with minimum rainfall in each 120kmx120km box|
| G | Afternoon_datetime_yyyymm_24km120km_iterate_filterRX1.txt | read in date and time when maximum rainfall occurs |

#### Output files:
|  | File names | Description |
| ----- | ------ | ------ |
| A	| Obtain_SMpercentiles_AfternoonRain_mmm_5degree1d_filterRX1.nc | record the spatial and temporal preferences when maximum rainfall occurs, contributed by soil moisture sourced from preceeding MCS or non-MCS rainfall |

-------
### <a name="head4"></a>4. ObtainSH_SMpercentile_NightRain_5degree1d.py

The script is similar to "ObtainSH_SMpercentile_AfternoonRain_5degree1d.py" but calculates the spatial and temporal preferences when nighttime rainfall occurs.

-------
### <a name="head5"></a>5. Obtain_WarmSM5ET5_correlations_160GridSmooth_detrend.py

The script calculates the coupling strength between initial soil moisture and its following pentad evapotranspiration, sourced from preceeding MCS or non-MCS rainfall. Note that the data has been smoothed and detrended.

#### Input files:
|  | File names | Description |
| ----- | ------ | ------ |
| A	| Obtain_ET_allpentads_yyyy_160GridSmooth_corrected.nc | read in pentad evapotranspiration of each year, with the components contributed by earlier-season MCS or non-MCS rainfall seperated |
| B | Obtain_SM123_allpentads_yyyy_160GridSmooth.nc | read in soil moisture before each pentad, with the components contributed by earlier-season MCS or non-MCS rainfall seperated |

#### Output files:
|  | File names | Description |
| ----- | ------ | ------ |
| A	| Obtain_mmmSM5_ET5_correlations_160GridSmooth_detrend.nc | record the coupling strength between initial soil moisture and following pentad evapotranspiration, contributed by MCS or non-MCS rainfall |

--------
### <a name="head6"></a>6. Obtain_WarmSM5PP5_correlations_160GridSmooth_detrend.py

The script calculates the coupling strength between initial soil moisture and a pentad (5-day) lagged pentad rainfall, with MCS or non-MCS components. Note that the data has been smoothed and detrended.

#### Input files:
|  | File names | Description |
| ----- | ------ | ------ |
| A	| Obtain_PP_allpentads_yyyy_160GridSmooth.nc | read in pentad rainfall of each year, with the components of MCS or non-MCS rainfall seperated |
| B | Obtain_SM123_allpentads_yyyy_160GridSmooth.nc | read in soil moisture before each pentad, with the components contributed by earlier-season MCS or non-MCS rainfall seperated |

#### Output files:
|  | File names | Description |
| ----- | ------ | ------ |
| A	| Obtain_mmmSM5_PP5_correlations_160GridSmooth_detrend.nc | record the coupling strength between initial soil moisture and a-pentad-lagged pentad rainfall, with MCS or non-MCS components |

--------
### <a name="head7"></a>7. Obtain_WarmSM5ET5xET5PP5_correlations_160GridSmooth_detrend.py

The script calculates the coupling strength between initial soil moisture and a pentad (5-day) lagged pentad rainfall, calculated by the production of coupling strength between soil moisture and following pentad evapotranspiration and the coupling strength between pentad evapotranspiration and a-pentad-lagged pentad rainfall. Each of them has MCS or non-MCS components. Note that the data has been smoothed and detrended.

#### Input files:
|  | File names | Description |
| ----- | ------ | ------ |
| A	| Obtain_ET_allpentads_yyyy_160GridSmooth_corrected.nc | read in pentad evapotranspiration of each year, with the components contributed by earlier-season MCS or non-MCS rainfall seperated |
| B | Obtain_SM123_allpentads_yyyy_160GridSmooth.nc | read in soil moisture before each pentad, with the components contributed by earlier-season MCS or non-MCS rainfall seperated |
| C	| Obtain_PP_allpentads_yyyy_160GridSmooth.nc | read in pentad rainfall of each year, with the components of MCS or non-MCS rainfall seperated |

#### Output files:
|  | File names | Description |
| ----- | ------ | ------ |
| A	| Obtain_mmmSM5ET5XET5PP5_correlations_160GridSmooth_detrend.nc | record another way to calculate coupling strength between initial soil moisture and a-pentad-lagged pentad rainfall, through the production of the two-legged coupling strengths |

-------
### <a name="head7"></a>Contact us
Huancui Hu: Huancui.Hu@pnnl.gov




