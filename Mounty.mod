###############################################################
#PCraster evapotranspiration model designed for use in montane environments
#Oliver Perkins July 2018

###############################################################

#! --radians

binding

### map inputs
 
 LULC = LULC.map; #model works with LC classification as per the ESA 20m land cover map http://2016africalandcover20m.esrin.esa.int/
 Dem = DEM.map; #m
 clone = Bool.map;
  
### time series inputs
 sunlight_hours_series = Sunlight_hours.tss;
 mintempseries = Min_temp.tss; #Celsuius
 maxtempseries = Max_temp.tss; #Celsius
 windseries = Wind.tss; #mph
 Rainseries = Rain.tss; #mm

timer

 1 365 1;
 reportdefault = 95+180..endtime;

initial

 PI = 3.1415926;
 e = 2.7182818;
  

 ### spatio-temporal constants
 day = 1; #(Julian Calendar)
 Latitude = -0.05; # degrees
 Daylight_hours = 12; #hours, see Allen et al. (1998) annex 2 
 Altitude_of_temperature_measurements = 915; #metres ASL
 
 
 
###################################################################
# Lookup tables
###################################################################

 ###Land cover properties
 
 G = lookupscalar(G.tbl, LULC); 
 Albedo = lookupscalar(Albedo.tbl, LULC); 
 Height = lookupscalar(Height.tbl, LULC); 
 LAI = lookupscalar(LAI.tbl, LULC); 
 
 ### Soil Hydrology
 
 ### Saturated, residual and normalised soil water content
 Theta_s = lookupscalar(Theta_s.tbl, LULC); #mm3/mm3
 Theta_r = lookupscalar(Theta_r.tbl, LULC); #mm3/mm3
 
 ### Soil Saturated Hydraulic Conductivity
 K_sat = lookupscalar(Ks.tbl, LULC); #mm day-1 (Wang et al., 2009)
  
 ### Soil water release curve parameters
 Alpha = lookupscalar(Van_genuchten_a.tbl, LULC); #Dimensionless
 N = lookupscalar(Van_genuchten_n.tbl, LULC); #Dimensionless
 M = 1 - 1/N; #Dimensionless 

 
 
###################################################################
# Constants for net energy balance
###################################################################

 ### Extra Terrestrial Solar radiation – (Allen et al., 1998)
 Latitude_Radians = (PI/180)*Latitude; #in Radians
 Solar_Constant = 0.0820; #MJ m-2 min-1
 Latent_Heat = 2.45; #MJ/kg
 
 ### Ground level radiation - see Paulescu et al., 2016
 Angstrom_Prescott_cloudy = 0.19840; #empirical coefficient 
 Angstrom_Prescott_clear  = 0.55161; # empirical coefficient
 Angstrom_Prescott_altitude = 0.04571*(Dem/1000); # impact of topography

 ### Net longwave radiation 
 Stefan_Boltzmann = 4.903 * (10**-9); #MJ K-4 m-2 day-1

 
#################################################################### 
# Constants for atmospheric pressure
###################################################################
 
 
 ### Atmospheric constants 
 
 ### Atmospheric Pressure
 P = 101.3*((293-(0.0065*Dem))/293)**5.26; #kPa
 
 ### Psychrometric constant 
 Gamma = (0.665*(10**-3))*P; #Dimensionless
 
 ### Specific Heat Capacity of Air 
 Cp = 1.004; #KJ kg-1 see - https://www.ohio.edu/mechanical/thermo/property_tables/air/air_cp_cv.html
 
 ### Correct Temperature for altitude
 Altitude_Temperature_Correction = (((Altitude_of_temperature_measurements - Dem )/ 304.8) *  1.98); # Change in Degrees, ICAO 1993

 
###################################################################
# Constants for Canopy structure
###################################################################


 ###Fractional Vegetation Cover as function of LAI (Beer's law; Jetten 2005)
 Light_Extinction_Coefficent = 0.4; #Dimensionless
 FVC = 1-(exp(-Light_Extinction_Coefficent*LAI)); 
  
 ### Canopy interception capacity as function of LAI (Vegas Galdos et al. 2012)
 Canopy_max = 2*ln(1+LAI); #mm
 Canopy_store = 0.7*Canopy_max; #mm 

 
###################################################################
# Constants for Soil Hydrology
################################################################### 
 
 
 ### Plant water stress parameters (Jetten 2005; Van Genuchten 1980)
 H50 = 8;  #Metres pressure head
 Soil_p = 3; #Dimensionless

 ### Initial Water content
 Theta = 0.9*Theta_s; # Dimensionless
 Theta_normalised = (Theta - Theta_r) / (Theta_s - Theta_r); #Dimensionless
 
 ###Total and initial soil storage capacity (Speich et al., 2018)
 Soil_Storage_Capacity = 2000 * Theta_s; #mm
 Soil_water_storage = Theta * 2000; #mm

 
 
dynamic

###################################################################
# Meteorological inputs
###################################################################


 ###Temperature corrected for altitude 
 min_temperature = timeinputscalar(mintempseries, 1) + Altitude_Temperature_Correction; #all in Celsius
 max_temperature = timeinputscalar(maxtempseries, 1) + Altitude_Temperature_Correction;
 mean_temperature = (max_temperature+min_temperature)/2;
 
 Wind = timeinputscalar(windseries, 1); #m/s
 Rain = timeinputscalar(Rainseries, 1); #mm 
 sunlight_hours = timeinputscalar(sunlight_hours_series, 1); #hours/day 
 
 max_temp_Kelvin = max_temperature + 273.15; # Kelvin needed for lw radiation 
 min_temp_Kelvin = min_temperature + 273.15; # Kelvin needed for lw radiation 
 mean_temp_Kelvin = mean_temperature + 273.15; # Kelvin for Air Density 
 
 
###################################################################
# Daily net energy balance
###################################################################


 ### Daily extraterrestrial radiation
 
 Decimation = 0.409*sin((2*PI)/365*day-1.39); #radians
 Earth_Sun_distance = cos((2*PI)*day/365)*0.033 + 1;  
 Sunset_angle = scalar(acos(-tan(Latitude_Radians) * tan(Decimation))); #radians
 Daily_rad_extra = ((24*60)/PI)*Solar_Constant*Earth_Sun_distance*(Sunset_angle*sin(Latitude_Radians)*sin(Decimation)+cos(Latitude_Radians)*cos(Decimation)*sin(Sunset_angle)); #MJ m-2 day-1
  
 ### Daily ground level radiation (Paulescu et al. 2016)
 Ground_radiation = (Angstrom_Prescott_cloudy + Angstrom_Prescott_clear*(sunlight_hours/Daylight_hours) + Angstrom_Prescott_altitude)*Daily_rad_extra; # MJ m-2 day-1

 ### albedo 
 Albedo_daily = Albedo*Ground_radiation; ##MJ m-2 day-1
  

 ### long wave radiation
 
 Actual_vapour_pressure = 0.6108*(e**((17.27*min_temperature)/(min_temperature + 237.3))); #KPa, assumes Tdew ~ Tmin
 Clearsky_radiation = (Angstrom_Prescott_cloudy + Angstrom_Prescott_clear + Angstrom_Prescott_altitude)*Daily_rad_extra; #MJ m-2 day-1
 Cloudfactor = Ground_radiation/Clearsky_radiation; #Dimensionless
 longwave = (((Stefan_Boltzmann*max_temp_Kelvin**4) + (Stefan_Boltzmann*min_temp_Kelvin**4))/2)*(0.34-0.14*sqrt(Actual_vapour_pressure))*(1.35*Cloudfactor - 0.35); #MJ m-2 day-1
  
 ### net shortwave radiation (vegetation)
 Net_Radiation = (Ground_radiation - Albedo_daily - longwave); #MJ m-2 day-1
 Rn_vegetation = Net_Radiation*FVC; ##MJ m-2 day-1
 
 ### Daily G value 
 Rn_subtractG = Net_Radiation-(G*Net_Radiation); ##MJ m-2 day-1
  

###################################################################
# Daily Atmospheric Variables (Allen et al. 1998)
###################################################################
 
 
 ### Saturated vapour pressure
 
 Saturatedpressure_mintemp = 0.6108*(e**((17.27*min_temperature)/(min_temperature + 237.3))); #kPa
 Saturatedpressure_maxtemp = 0.6108*(e**((17.27*max_temperature)/(max_temperature + 237.3))); #kPa
 Saturatedpressure_meantemp = 0.6108*(e**((17.27*mean_temperature)/(mean_temperature + 237.3))); #kPa
 Saturated_Pressure = (Saturatedpressure_maxtemp+Saturatedpressure_mintemp)/2; #kPa
 
 ### Vapour_Pressure_Deficit
 Vapour_Pressure_Deficit = Saturated_Pressure - Actual_vapour_pressure; #kPa

 ### Derivative of Saturation pressure for Daily mean temperature
 Delta = (4098*Saturatedpressure_meantemp)/((mean_temperature+237.3)**2); #d kPa by d Celsius

 ### Density of Dry Air
 Rho = (1000*P) / (287.058 * mean_temp_Kelvin) ; #kg/m3; Ideal Gas Equation; also converts kPa to Pascals

 ### Aerodynamic resistance - equation from FAO
 ### Set measurement height to 5 metres to reflect high canopy
  
 R_aero = (ln((5-0.67*Height)/(0.123*Height)) * ln((5-0.67*Height)/(0.1*0.123*Height))) / ((0.41**2)*(Wind/0.838)); #parameterised as per Allen et al. 1998, windconversion factor from Table Annex 2
 Raero_soil = (ln((5-0.67*0.001)/(0.123*0.001)) * ln((5-0.67*0.001)/(0.1*0.123*0.001))) / ((0.41**2)*(Wind/0.838)); #for soil evaporation, parameters from Mohammed et al., 1997

 
###################################################################
# Hydrological cycle
###################################################################
  
  
 ###Interception
 Canopy_tminus1 = Canopy_store; #for calculating Ground Rain
 Canopy_store = min(Canopy_max, Canopy_store + ((Canopy_max - Canopy_store)*(1-exp(-1/Canopy_max*Rain)))); #mm
 
 ### Wet fraction of canopy (after Speich et al. 2018)
 Wet_fraction = (Canopy_store / Canopy_max) ** 2/3; #Dimensionless
 
 ### Ground rain
 Ground_rain = Rain - (Canopy_store - Canopy_tminus1); #mm
 
 
###################################################################
# Soil Hydrology
###################################################################


 ### update Soil water balance 
 Soil_water_tminus1 = Soil_water_storage; #for calculating runoff
 Soil_water_storage = min(Soil_Storage_Capacity, Soil_water_storage + Ground_rain); #mm
 runoff = Ground_rain - (Soil_water_storage - Soil_water_tminus1); #mm, for verification
 
 ### Update Theta
 Theta = Soil_water_storage / 2000; #mm3/mm-3
 Theta_normalised = (Theta - Theta_r) / (Theta_s - Theta_r); #Dimensionless
 
 ### Percolation 
 Percolation = if(Theta_normalised gt 0.9, K_sat*((Theta - Theta_r)**3.5), 0); #mm - (after Mualem et al., 1976)
  
 ### Matric potential (Van Genuchten 1980)
 H = (((1/Theta_normalised)**1/M - 1)**(1/N)) / Alpha; #Metres
 
 ### Vegetation water stress function  (Van Genuchten 1987; cited in Yang et al. 2013)
 Water_stress_function = (1 / (1 + (H / H50)**Soil_p)); #Dimensionless
 
 ### Minimum Canopy Resistance (Allen et al., 1998)
 R_canopy = 100/((LAI/2) * Height); #s m-1 
 
 ### Soil evaporative resistance (Kang et al., 2009)
 Rc_soil = exp(8.2-(4.225*Theta)/Theta_s); #s m-1

 
###################################################################
# Calculate evapotranspiration
###################################################################


 ### Penman-Monteith potential evapotranspiration (Monteith 1965)
 LE = ((Delta*Net_Radiation) + ((Rho*Cp*Vapour_Pressure_Deficit) /R_aero))/(Delta + (1+R_canopy/R_aero)*Gamma); #MJ m-2 day-1
 ETp = LE / Latent_Heat; #conversion to mm
 report ETpotential.tss = timeoutput(clone, ETp); #mm
 
### Interception Evaporation (after Speich et al., 2018)
 Rn_Int = Rn_vegetation * Wet_fraction; #MJ m-2 day-1
 E_Int_potential = ((Delta*Rn_Int + (FVC * (Rho * Cp * Vapour_Pressure_Deficit)/R_aero)) / (Delta + Gamma)); #MJ m-2 day-1
 E_Int_potential = E_Int_potential / Latent_Heat; #mm
 report E_Int = min(Canopy_store, E_Int_potential); #mm
 report E_Int.tss = timeoutput(clone, E_Int); #mm
 
 ### Transpiration (simplified from Speich et al: uses wet/dry fraction not energy balance)
 Rn_trans = Rn_vegetation * (1-Wet_fraction); ## #MJ m-2 day-1
 LE_trans = ((Delta*Rn_trans + (FVC * (Rho * Cp * Vapour_Pressure_Deficit)/R_aero)) /(Delta + (1+R_canopy/R_aero)*Gamma)) * Water_stress_function; #MJ m-2 day-1
 report Transpiration = LE_trans / Latent_Heat; #mm -d 
 report Transpiration.tss = timeoutput(clone, Transpiration); #mm
 
 ###Soil Evaporation
 Rn_soil = Rn_subtractG*(1-FVC); ##MJ m-2 day-1
 E_Soil = ((Delta*Rn_soil + ((1-FVC)*(Rho * Cp * Vapour_Pressure_Deficit)/ Raero_soil)) / (Delta + (1+Rc_soil / Raero_soil)*Gamma)) * Theta_normalised; #MJ m-2 day-1
 report E_Soil = E_Soil / Latent_Heat; #mm
 report E_Soil.tss = timeoutput(clone, E_Soil); #mm
 
 ### Total actual ET
 report Total_ET = E_Int + Transpiration + E_Soil; #mm
 report ET_Total.tss = timeoutput(clone, Total_ET); #mm
 
 ### Average by land use class
 report forest_total.tss = timeoutput(LULC eq 1, Total_ET); #mm 
 report agriculture_total.tss = timeoutput(LULC eq 4, Total_ET); #mm



#################################################################### 
# Variable updates
###################################################################
 
 
 ### update water balance
 Canopy_store = max(0, Canopy_store - E_Int); #mm
 Soil_water_storage = max(Theta_r*2000, Soil_water_storage - Transpiration - E_Soil - Percolation); #mm
 
 ### Update day
 day = day + 1;
  
  
###################################################################
# Check Energy & Water balance
###################################################################

 ### Energy (Wang and Dickinson 2012)
report Energy_balance.tss = timeoutput(clone, Ground_radiation - (Rn_trans + Rn_Int + Rn_soil + (Net_Radiation*G*(1-FVC)) + longwave + (Albedo*Ground_radiation))); #MJ m-2 day-1

 ### Water (Wang and Dickinson 2012)
 report Water_balance.tss = timeoutput(clone, Rain - (Transpiration + E_Soil + E_Int + Percolation + runoff + (Canopy_store - Canopy_tminus1) + (Soil_water_storage - Soil_water_tminus1))); #mm


