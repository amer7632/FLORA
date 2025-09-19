###########################################################################
####    pyFLORA - pythonic FAST LAND OCCUPANCY AND REACTION ALGORITHM   ###
####                       Paleoclimate runs                            ###
###########################################################################


### SETUP ###

import numpy as np
import scipy.io as sp
from scipy.stats import norm

file_to_open =f'./data/paleo2021_data1.mat'
flora_interpstack = sp.loadmat(file_to_open)
flora_lat = flora_interpstack['lat']
flora_lon = flora_interpstack['lon']
#contains
#CO2_data = CO2 (ppm) value
            #(0) 10, (1) 50, (2) 100, (3) 160, (4) 180, (5) 200, (6) 240,
            #(7) 280, (8) 360, (9) 420, (10) 560, (11) 700, (12) 840, 
            #(13) 1120, (14) 1400, (15) 1680, (16) 2240, (17) 2800, 
            #(18) 3360, (19) 4200, (20) 5600, (21) 7000, (22) 14000, 
            #(23) 28000
#land_data = true/false where land is present
#lat/lon = lat and lon coordinates
#runoff_data = runoff for each time point and CO2 concentration
#temp_data = temperature for each time point and CO2 concentration
#time_data = time (Ma) values
            #(1) 540, (2) 500, (3) 470, (4) 450, (5) 430, (6) 400, (7) 370,
            #(8) 340, (9) 300, (10) 280, (11) 260, (12) 245, (13) 220, 
            #(14) 200, (15) 180, (16) 145, (17) 90, (18) 70, (19) 52, 
            #(20) 30, (21) 15, (22) 0

# CO2 and O2 levels
co2_level = np.asarray([21, 19, 19, 18, 17, 16, 15, 12, 11, 12, 12, 14, 16, 16, 16, 13, 14, 11, 12, 11, 9, 8])
co2_m_stdev = np.asarray([20, 18, 18, 17, 16, 12, 10, 6, 3, 8, 7, 10, 11, 12, 11, 10, 10, 4, 10, 6, 4, 4])
co2_p_stdev = np.asarray([22, 20, 20, 19, 18, 19, 17, 14, 13, 14, 14, 16, 18, 18, 17, 16, 17, 14, 14, 13, 11, 11]) 
o2 = np.asarray([ 1.9, 5.14, 3.74, 5.05, 14.94, 16.76, 16.22, 24.76, 25.52, 27.25, 31.05, 30.77, 32.38, 28.6, 25.13, 25.86, 26.65, 27.69, 25.12, 21.41, 21.57, 20.72])


# Number of loops
nsteps = 200
# Longitude and latitude
x_lon = 40
y_lat = 48

# starting time
t = 0 #1 in matlabe
# timesteps in 1 year
timestep = 1 

# reservoirs
# Reservoirs #
C_leaf_tem = np.full((x_lon, y_lat, nsteps), np.nan)
C_leaf_bor = np.full((x_lon, y_lat, nsteps), np.nan)
C_leaf_tro = np.full((x_lon, y_lat, nsteps), np.nan)
R_leaf_tem = np.full((x_lon, y_lat, nsteps), np.nan)
R_leaf_bor = np.full((x_lon, y_lat, nsteps), np.nan)
R_leaf_tro = np.full((x_lon, y_lat, nsteps), np.nan)
leaf_turnover_tem = np.full((x_lon, y_lat, nsteps), np.nan)
leaf_turnover_bor = np.full((x_lon, y_lat, nsteps), np.nan)
leaf_turnover_tro = np.full((x_lon, y_lat, nsteps), np.nan)
NPP_tem = np.full((x_lon, y_lat, nsteps), np.nan)
NPP_bor = np.full((x_lon, y_lat, nsteps), np.nan)
NPP_tro = np.full((x_lon, y_lat, nsteps), np.nan)
biomass_tem = np.full((x_lon, y_lat, nsteps), np.nan)
biomass_bor = np.full((x_lon, y_lat, nsteps), np.nan)
biomass_tro = np.full((x_lon, y_lat, nsteps), np.nan)

### some parameters ###
# Q10 temperature dependence
t25 = 2600
q_10t = 0.57

# Intracellular partial pressure 
v = 0.7
kc25 = 30
ko25 = 3e4
q_10c = 2.1
q_10o = 1.2
s = 0.015
alpha = 0.08
theta = 0.7
lr_max = 0.75

# Leaf c:n ratio
CN_leaf = 29

# Tissue respiration rate at 10 degree C
r_tem = 0.055 * 365 # 0.066 for some PFTs; gC/gN/d -> gC/gN/year
r_bor = 0.066 * 365
r_tro = 0.011 * 365

# Growth respiraition
R_growth = 0.25

# Leaf longevity; ranges from 0.5 - 1 depending on type of plant
life_leaf_tem = 0.75 
life_leaf_bor = 0.75
life_leaf_tro = 1


#Insolation 
ins_present = 150 + 250 * np.outer(norm.pdf(flora_lat, 0, 40) / norm.pdf(0, 0, 40), np.ones(x_lon))

for flora_time_ind, flora_time in enumerate(flora_times[:1]):
    print(flora_time_ind, flora_time)

    # we will count forward in time based on our climate data
    # get index of co2
    co2_ind = co2_level[flora_time_ind]
    co2_ppm = flora_co2_data[co2_ind] # ppm = micromol

    # land
    land = flora_land_data[:, :, flora_time_ind]
    land[land == 0] = np.nan

    # temperature, masked by land
    temp = flora_temp_data[:,:, co2_ind, flora_time_ind] * land

    # runoff, masked by land
    runoff = flora_runoff_data[:,:, co2_ind, flora_time_ind] * land

    # mask out where ice forms?
    temp[temp < -10] = np.nan
    

    #Ambient partial pressure
    po2 = o2[flora_time_ind] * 1000
    pco2 = co2_ppm

    tf = t25*(q_10t**((temp - 25)*0.1))
    tstar = po2/(2*tf )
    pi = v * pco2
    kc = kc25*(q_10c**((temp - 25)*0.1))
    ko = ko25*(q_10o**((temp - 25)*0.1))
    c2 = (pi - tstar)/(pi + kc**(1 + (po2/ko)))

    # make ftemps for temperate, boreal and tropical regions
    ftemp_tem = norm.pdf(temp, loc=15, scale=20)
    ftemp_bor = norm.pdf(temp, loc=0, scale=20)
    ftemp_tro = norm.pdf(temp, loc=27, scale=7)
    c1_bor = alpha*ftemp_bor**((pi - tstar)/(pi + 2*tstar)) 
    c1_tem = alpha*ftemp_tem**((pi - tstar)/(pi + 2*tstar)) 
    c1_tro = alpha*ftemp_tro**((pi - tstar)/(pi + 2*tstar))

    sigma = (1 - ((c2 - s)/(c2 - theta*s)))**0.5


    here
    
        #Change in insolation overtime
        ins( : , : ) = ins_present - (ins_present * 4.6/100 * (abs(time_data(a))/570)) ; 
        
        # Maintenance respiration # 

        #Arrhenius equation; temperature function for respiration
        g_T ( : , : ) = exp( 308.56 .* ( ( 1 / 56.02) - ( 1 ./ ( tmp_avg( : , : ) + 46.02 ) ) ) ) ;

        ### Water presence ###

        #converting actual runoff into a sigmoidal curve; water presence ranked
        #from 0-1
        #test_water = runoff_re ; 
        #test_water(test_water == 0 | isnan(test_water)) = NaN ; 
        #test_water = circshift( test_water, -6 ) ; 
        water_stress = 1 - (1 ./ ( 1 + exp(0.005 .* (runoff_re - 450)))) ;

        ### PHOTOSYNTHESIS ###
        photosynth_tem( : , : ) = 10 * 365 * ins( : , : ) .* ( c1_tem( : , : ) ./ c2( : , : ) ) .* ( c2( : , : ) - ( 2 * theta - 1 ) * s - 2 .* ( c2( : , : ) - theta * s ) .* sigma( : , : ) ) .* water_stress( : , : ) ; 
        photosynth_bor( : , : ) = 10 * 365 * ins( : , : ) .* ( c1_bor( : , : ) ./ c2( : , : ) ) .* ( c2( : , : ) - ( 2 * theta - 1 ) * s - 2 .* ( c2( : , : ) - theta * s ) .* sigma( : , : ) ) .* water_stress( : , : ) ; 
        photosynth_tro( : , : ) = 10 * 365 * ins( : , : ) .* ( c1_tro( : , : ) ./ c2( : , : ) ) .* ( c2( : , : ) - ( 2 * theta - 1 ) * s - 2 .* ( c2( : , : ) - theta * s ) .* sigma( : , : ) ) .* water_stress( : , : ) ; 
        # gC/ m^2/ year

        ### CARBON IN LEAF ###
        C_leaf_tem( : , : , 1 ) = lr_max .* photosynth_tem( : , : ) ;
        C_leaf_bor( : , : , 1 ) = lr_max .* photosynth_bor( : , : ) ;
        C_leaf_tro( : , : , 1 ) = lr_max .* photosynth_tro( : , : ) ;
        # gC/ m^2/ year

        ### BIOMASS ###
        biomass_tem( : , : , 1 ) = 2.5e4 .* land( : , : ) ; 
        biomass_tem( biomass_tem == 0 ) = NaN ; 
        biomass_bor( : , : , 1 ) = 2.5e4 .* land( : , : ) ; 
        biomass_bor( biomass_bor == 0 ) = NaN ; 
        biomass_tro( : , : , 1 ) = 2.5e4 .* land( : , : ) ; 
        biomass_tro( biomass_tro == 0 ) = NaN ; 
        # gC/m^2/year

       #Biomass Calculation#
       for n = 1 : nsteps
    
            ### Leaf respiration ###
            R_leaf_bor( :, :, n ) = r_bor * ( C_leaf_bor( :, :, n ) / CN_leaf ) .* g_T( : , : ) ;
            R_leaf_bor(R_leaf_bor < 0) = 0 ; 
            R_leaf_tem ( :, :, n ) = r_tem * ( C_leaf_tem( :, :, n ) / CN_leaf ) .* g_T( : , : ) ;
            R_leaf_tem(R_leaf_tem < 0) = 0 ; 
            R_leaf_tro ( :, :, n ) = r_tro * ( C_leaf_tro( :, :, n ) / CN_leaf ) .* g_T( : , : ) ;
            R_leaf_tro(R_leaf_tro < 0) = 0 ; 
        
            ### NPP ##
            NPP_bor ( :, :, n ) = ( 1 - R_growth ) .* ( photosynth_bor( : , : ) - R_leaf_bor( :, :, n ) ) ;
            NPP_tem ( :, :, n ) = ( 1 - R_growth ) .* ( photosynth_tem( :, : ) - R_leaf_tem( :, :, n ) ) ; 
            NPP_tro ( :, :, n ) = ( 1 - R_growth ) .* ( photosynth_tro( : , : ) - R_leaf_tro( :, :, n ) ) ; 
            #gC/m2/year
        
            ### Carbon in leaf allocation ###
            C_leaf_bor( :, :, n+1 ) = ( C_leaf_bor( :, :, n ) .* ( 1 - life_leaf_bor ) ) + ( lr_max .* NPP_bor( :, :, n ) )  ;
            C_leaf_tem( :,:, n+1 ) = ( C_leaf_tem( :, :, n ) .* ( 1 - life_leaf_tem ) ) + ( lr_max .* NPP_tem( :, :, n ) )  ;
            C_leaf_tro( :, :, n+1 ) = ( C_leaf_tro( :, :, n ) .* ( 1 - life_leaf_tro ) ) + ( lr_max .* NPP_tro( :, :, n ) )  ;
        
            ### Biomass ###
            biomass_bor( :, :, n+1 ) = biomass_bor( :, :, n ) + ( C_leaf_bor( :, :, n ) - 0.1 * biomass_bor( :, : , n ) ) .* time_step ; 
            biomass_tem( :, :, n+1 ) = biomass_tem( :, :, n ) + ( C_leaf_tem( :, :, n ) - 0.1 * biomass_tem( :, :, n ) ) .* time_step ; 
            biomass_tro( :, :, n+1 ) = biomass_tro( :, :, n ) + ( C_leaf_tro( :, :, n ) - 0.1 * biomass_tro( :, : , n ) ) .* time_step ; 

       end

        final_biomass_tem( :, :, a ) = biomass_tem( :, :, end ) ; 
        final_biomass_bor( :, :, a ) = biomass_bor( :, :, end ) ; 
        final_biomass_tro( :, :, a ) = biomass_tro( :, :, end ) ;
        
        final_NPP_tem( :, :, a ) = NPP_tem( :, :, end ); 
        final_NPP_bor( :, :, a ) = NPP_bor( :, :, end ); 
        final_NPP_tro( :, :, a ) = NPP_tro( :, :, end );
        
        #Saving temp info 
        temp_end( :, :, a ) = tmp_avg ; 
        runoff_end( :, :, a ) = runoff_re ; 
end



###########################
###     Competition     ###
###########################

for a = 1 : 22
    for i = 1 : x_lon
        for j = 1 : y_lat
           final_biomass( i, j, a ) =  max( final_biomass_tem( i, j, a ),max( final_biomass_bor( i, j, a ), final_biomass_tro( i, j, a ) ) ) ;
           if final_biomass( i, j, a ) == final_biomass_tem( i, j, a )
               biome( i, j, a ) = 1 ; 
               final_NPP(i,j,a) = final_NPP_tem(i,j,a) ;
           elseif final_biomass( i, j, a ) == final_biomass_bor( i, j, a )
               biome( i, j, a ) = 2 ; 
               final_NPP(i,j,a) = final_NPP_bor(i,j,a) ; 
           elseif final_biomass( i, j, a ) == final_biomass_tro( i, j, a )
               biome( i, j, a ) = 3 ;
               final_NPP(i,j,a) = final_NPP_tro(i,j,a) ; 
           else
               biome( i, j, a ) = 4 ; 
               final_NPP(i,j,a) = NaN ; 
           end
        end
    end
end



###########################
###    Save results     ###
###########################

#collects necessary end result data
# save('FLORA_Paleorun', 'CO2_m_stdev', 'CO2_p_stdev', 'final_biomass', 'biome', 'temp_end', 'runoff_end', 'final_biomass_bor', 'final_biomass_tem', 'final_biomass_tro', 'CO2_level', 'O2', 'time', 'x_lon', 'y_lat', 'final_NPP' )



