%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    FLORA - FAST LAND OCCUPANCY AND REACTION ALGORITHM   %%%
%%%                 Present day validation                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% SETUP %%%

load FLORA_present_day_data.mat
%contains
%CDIAC_biomass: IPCC biomss with reduced resolution to match proxy data
%land_cru, lat_cru, lon_cru: CRU data
%runoff: runoff from Fekete et al. 2000
%tmp_avg: CRU temperatue data


%%% Arrays %%%

% Number of loops %
nsteps = 157 ; 
% Longitude and latitude %
x_lon = 360 ; 
y_lat = 720 ; 

% Reservoirs %
C_leaf_tem = NaN( x_lon, y_lat, nsteps ) ;
C_leaf_bor = NaN( x_lon, y_lat, nsteps ) ; 
C_leaf_tro = NaN (x_lon, y_lat, nsteps ) ; 
R_leaf_tem = NaN( x_lon, y_lat, nsteps ) ;
R_leaf_bor = NaN( x_lon, y_lat, nsteps ) ;
R_leaf_tro = NaN( x_lon, y_lat, nsteps ) ;
leaf_turnover_tem = NaN( x_lon, y_lat, nsteps ) ;
leaf_turnover_bor = NaN( x_lon, y_lat, nsteps ) ;
leaf_turnover_tro = NaN( x_lon, y_lat, nsteps ) ;
NPP_tem = NaN( x_lon, y_lat, nsteps ) ;
NPP_bor = NaN( x_lon, y_lat, nsteps ) ;
NPP_tro = NaN( x_lon, y_lat, nsteps ) ;
biomass_tem = NaN( x_lon, y_lat, nsteps ) ;
biomass_bor = NaN( x_lon, y_lat, nsteps ) ;
biomass_tro = NaN( x_lon, y_lat, nsteps ) ;

% Ice = no land %
for i = 1 : x_lon
    for j = 1 : y_lat   
        if tmp_avg( i, j ) < -10
            tmp_avg( i, j ) = NaN ; 
        else
        end
    end
end

t( 1 ) = 1 ;
time_step = 1 ; %year


%%% CONSTANTS %%%

%Ambient partial pressure
pO2 = 20.9 * 1000  ; %Pa
pCO2 =  280 ; %Pa %340 used in Sitch et al. 2003

%Q10 temperature dependence
t25 = 2600 ; 
q_10t = 0.57 ; 
tf( : , : ) = t25 * ( q_10t .^ ( ( tmp_avg( : , : ) - 25 ) * 0.1 ) ) ;
tstar( : , : ) = pO2 ./ ( 2 * tf ) ; 
% tstar = CO2 compensation point 
% pO2 = ambient partial pressure of O2
% tf = specificity factor that reflets the ability of Rubisco to discriminate between CO2 and O2 
% t_25 = tf value at 25oC 
% q_10t = temperature sensitivity parameter 

%Intracellular partial pressure 
v = 0.7 ; 
% P = 100 * 1000 ; 
% ca = 355 ; %umol/mol
% pi = v * P * ca ;
pi = v * pCO2 ; 
% pi = intracellular partial pressure of CO2; Pa
% P = atmospheric pressure; Pa
% ca = ambient CO2 concentration
% v = parameter with a positive value <= 0.8 (C3 plants); max value used
kc25 = 30 ; % Pa
ko25 = 3e4 ; % Pa
q_10c = 2.1 ; 
q_10o = 1.2 ; 
kc( : , : ) = kc25 * ( q_10c .^ ( ( tmp_avg( : , : ) - 25 ) * 0.1 ) ) ;
ko( : , : ) = ko25 * ( q_10o .^ ( ( tmp_avg( : , : ) - 25 ) * 0.1 ) ) ;
c2( : , : ) = ( pi - tstar( : , : ) ) ./ ( pi + kc( : , : ) .* ( 1 + pO2 ./ ko( : , : ) ) ) ; 
% kc and ko = kinetic parameters with a Q10 dependence on temperature; Pa

s = 0.015 ; 
% s = constant leaf respiration as a fraction of Rubsico capacity (C3 plants)

alpha = 0.08 ; %0.053 for C4 plants
theta = 0.7 ; 
ftemp_tem( : , : ) = normpdf( tmp_avg( : , : ), 15, 20 ) ;
ftemp_bor( : , : ) = normpdf( tmp_avg( : , : ), 0, 20 ) ;
ftemp_tro( : , : ) = normpdf( tmp_avg( : , : ), 27, 7 ) ; 
c1_bor( : , : ) = alpha * ftemp_bor .* ( ( pi - tstar( : , : ) ) ./ ( pi + 2 * tstar( : , : ) ) ) ; 
c1_tem( : , : ) = alpha * ftemp_tem .* ( ( pi - tstar( : , : ) ) ./ ( pi + 2 * tstar( : , : ) ) ) ; 
c1_tro( : , : ) = alpha * ftemp_tro .* ( ( pi - tstar( : , : ) ) ./ ( pi + 2 * tstar( : , : ) ) ) ;
% alpha = effective ecosystem-level quantum efficiency
% theta = co-limitation shape parameter
% ftemp = temperature inhibition function limiting photosynthesis at low and high temps

sigma( : , : ) = ( 1 - ( c2( : , : ) - s ) ./ ( c2( : , : ) - theta * s ) ) .^ 0.5 ;
ins( : , : ) = 150 + 250 .* normpdf( lat_cru, 0, 40 ) ./ normpdf( 0, 0, 40 ) .* ones( x_lon, y_lat ) ;
%insolation in place of PAR

% Maintenance respiration % 
lr_max = 0.75 ; 
%lr_max = leaf to root max ratio; 1 in nonstressed conditions, 0.75 for herbaceous

%Leaf c:n ratio
CN_leaf = 29 ;

%Tissue respiration rate at 10 degree C for Tropicals
r_tem = 0.055 * 365 ; %0.066 for some PFTs; gC/gN/d -> gC/gN/year
r_tro = 0.011 * 365 ; 
r_bor = 0.066 * 365 ; 

%Growth respiraition
R_growth = 0.25 ; 

%Modified Arrheinus equation 
g_T( : , : ) = exp( 308.56 .* ( ( 1 / 56.02) - ( 1 ./ ( tmp_avg( : , : ) + 46.02 ) ) ) ) ;

%Leaf turnover; ranges from 0.5 - 1 depending on type of plant
life_leaf_tem = 0.75 ; 
life_leaf_bor = 0.75 ; 
life_leaf_tro = 1 ; 

%%% Water presence %%%

%converting actual runoff into a sigmoidal curve; water presence ranked
%from 0-1
for i = 1 : 360
    for j = 1 : 720
        if runoff( i , j ) == 0 || runoff( i , j ) == NaN
            test_water( i , j ) = NaN; 
        else
            test_water( i , j ) = runoff( i , j ) ; 
        end
    end
end
test_water = circshift( test_water, -6 ) ; 

water_stress = sigmf( test_water, [ 0.005 450 ] ) ; 


%%% PHOTOSYNTHESIS %%%
photosynth_tem( : , : ) = 10 * 365 * ins( : , : ) .* ( c1_tem( : , : ) ./ c2( : , : ) ) .* ( c2( : , : ) - ( 2 * theta - 1 ) * s - 2 .*( c2( : , : ) - theta * s ) .* sigma( : , : ) ) .* water_stress( : , : ) ; 
photosynth_bor( : , : ) = 10 * 365 * ins( : , : ) .* ( c1_bor( : , : ) ./ c2( : , : ) ) .* ( c2( : , : ) - ( 2 * theta - 1 ) * s - 2 .*( c2( : , : ) - theta * s ) .* sigma( : , : ) ) .* water_stress( : , : ) ; 
photosynth_tro( : , : ) = 10 * 365 * ins( : , : ) .* ( c1_tro( : , : ) ./ c2( : , : ) ) .* ( c2( : , : ) - ( 2 * theta - 1 ) * s - 2 .*( c2( : , : ) - theta * s ) .* sigma( : , : ) ) .* water_stress( : , : ) ; 
% gC/ m^2/ year


%%% CARBON IN LEAF %%%
C_leaf_tem( : , : , 1 ) = lr_max .* photosynth_tem( : , : ) ;
C_leaf_bor( : , : , 1 ) = lr_max .* photosynth_bor( : , : ) ;
C_leaf_tro( : , : , 1 ) = lr_max .* photosynth_tro( : , : ) ;
% gC/ m^2/ year


%%% BIOMASS %%%
biomass_tem( : , : , 1 ) = 0.1 .* land_cru( : , : ) ; 
biomass_tem( biomass_tem == 0 ) = NaN ; 
biomass_bor( : , : , 1 ) = 0.1 .* land_cru( : , : ) ; 
biomass_bor( biomass_bor == 0 ) = NaN ; 
biomass_tro( : , : , 1 ) = 0.1 .* land_cru( : , : ) ; 
biomass_tro( biomass_tro == 0 ) = NaN ; 
% gC/m^2/year



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Biomass of PTFs   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Boreal plants%
tic
for n = 1 : nsteps
    
    %%% Leaf respiration %%%
    R_leaf_bor ( :, :, n ) = r_bor * ( C_leaf_bor( :, :, n ) / CN_leaf ) .* g_T( : , : ) ;
    for i = 1 : x_lon
        for j = 1 : y_lat
            if R_leaf_bor ( i , j , n ) > 0
               %do nothing
            else
                R_leaf_bor( i , j , n ) = 0 ;
            end
        end
    end
    %%% NPP %%
    NPP_bor ( :, :, n ) = ( 1 - R_growth ) .* ( photosynth_bor( : , : ) - R_leaf_bor( :, :, n ) ) ; 
    %gC/m2/year

    %%% Carbon in leaf allocation %%%
    C_leaf_bor( :, :, n+1 ) = ( C_leaf_bor( :, :, n ) .* ( 1 - life_leaf_bor ) ) + ( lr_max .* NPP_bor( :, :, n ) )  ;

    %%% Biomass %%%
    biomass_bor( :, :, n+1 ) = biomass_bor( :, :, n ) + ( C_leaf_bor( :, :, n ) - 0.1 * biomass_bor( :, : , n ) ) .* time_step ; 
    
end
toc

%Temperate plants%
tic
for n = 1 : nsteps

    %%% Leaf respiration %%%
    R_leaf_tem ( :, :, n ) = r_tem * ( C_leaf_tem( :, :, n ) / CN_leaf ) .* g_T( : , : ) ;
    for i = 1 : x_lon
        for j = 1 : y_lat
            if R_leaf_tem( i , j , n ) > 0
               %do nothing
            else
                R_leaf_tem( i , j , n ) = 0 ;
            end
        end
    end
    %%% NPP %%
    NPP_tem ( :, :, n ) = ( 1 - R_growth ) .* ( photosynth_tem( :, : ) - R_leaf_tem( :, :, n ) ) ; 
    %gC/m2/year

    %%% Carbon in leaf allocation %%%
    C_leaf_tem( :,:, n+1 ) = ( C_leaf_tem( :, :, n ) .* ( 1 - life_leaf_tem ) ) + ( lr_max .* NPP_tem( :, :, n ) )  ;

    %%% Biomass %%%
    biomass_tem( :, :, n+1 ) = biomass_tem( :, :, n ) + ( C_leaf_tem( :, :, n ) - 0.1 * biomass_tem( :, :, n ) ) .* time_step ; 

end
toc

%Tropical plants%
tic
for n = 1 : nsteps

    %%% Leaf respiration %%%
    R_leaf_tro ( :, :, n ) = r_tro * ( C_leaf_tro( :, :, n ) / CN_leaf ) .* g_T( : , : ) ;
    for i = 1 : x_lon
        for j = 1 : y_lat
            if R_leaf_tro ( i , j , n ) > 0
               %do nothing
            else
                R_leaf_tro( i , j , n ) = 0 ;
            end
        end
    end
    %%% NPP %%
    NPP_tro ( :, :, n ) = ( 1 - R_growth ) .* ( photosynth_tro( : , : ) - R_leaf_tro( :, :, n ) ) ; 
    %gC/m2/year

    %%% Carbon in leaf allocation %%%
    C_leaf_tro( :, :, n+1 ) = ( C_leaf_tro( :, :, n ) .* ( 1 - life_leaf_tro ) ) + ( lr_max .* NPP_tro( :, :, n ) )  ;

    %%% Biomass %%%
    biomass_tro( :, :, n+1 ) = biomass_tro( :, :, n ) + ( C_leaf_tro( :, :, n ) - 0.1 * biomass_tro( :, : , n ) ) .* time_step ; 

end
toc

%Biomass
final_biomass_tem = biomass_tem( :, :, end ) ; 
final_biomass_bor = biomass_bor( :, :, end ) ; 
final_biomass_tro = biomass_tro( :, :, end );
% gC/m2/year 



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Competition     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : x_lon
    for j = 1 : y_lat
       final_biomass( i , j ) =  max( final_biomass_tem( i , j ), max( final_biomass_bor( i , j ), final_biomass_tro( i , j ) ) ) ;
       %1 = temperate, 2 = boreal, 3 = tropical, 4 = ice/desert
       if final_biomass( i , j ) == final_biomass_tem( i , j )
           biome( i , j ) = 1 ; 
       elseif final_biomass(i,j) == final_biomass_bor( i , j )
           biome( i , j ) = 2 ; 
       elseif final_biomass(i,j) == final_biomass_tro( i , j )
           biome( i , j ) = 3 ;
       else
           biome( i , j ) = 4 ; 
       end
    end
end

biome = biome .* land_cru ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Comparison      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%comparing CDIAC carbon values to model 
c_biomass = CDIAC_biomass / 10 ; 
% kg/ha -> gC/m^2
compare_carbon = final_biomass - c_biomass ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Save results     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%collects necessary end result data
save('FLORA_validation', 'final_biomass', 'biome', 'c_biomass', 'compare_carbon', 'final_biomass_bor', 'final_biomass_tem', 'final_biomass_tro', 'x_lon', 'y_lat' )

