%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    FLORA - FAST LAND OCCUPANCY AND REACTION ALGORITHM   %%%
%%%                 Paleoclimate runs                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% SETUP %%%

clear
load data/paleo2021_data1.mat
%contains
%CO2_data = CO2 (ppm) value
            %(1) 10, (2) 50, (3) 100, (4) 160, (5) 180, (6) 200, (7) 240,
            %(8) 280, (9) 360, (10) 420, (11) 560, (12) 700, (13) 840, 
            %(14) 1120, (15) 1400, (16) 1680, (17) 2240, (18) 2800, 
            %(19) 3360, (20) 4200, (21) 5600, (22) 7000, (23) 14000, 
            %(24) 28000
%land_data = true/false where land is present
%lat/lon = lat and lon coordinates
%runoff_data = runoff for each time point and CO2 concentration
%temp_data = temperature for each time point and CO2 concentration
%time_data = time (Ma) values
            %(1) 540, (2) 500, (3) 470, (4) 450, (5) 430, (6) 400, (7) 370,
            %(8) 340, (9) 300, (10) 280, (11) 260, (12) 245, (13) 220, 
            %(14) 200, (15) 180, (16) 145, (17) 90, (18) 70, (19) 52, 
            %(20) 30, (21) 15, (22) 0

%CO2 and O2 levels
CO2_level = [21, 19, 19, 18, 17, 16, 15, 12, 11, 12, 12, 14, 16, 16, 16, 13, 14, 11, 12, 11, 9, 8];
CO2_m_stdev = [20, 18, 18, 17, 16, 12, 10, 6, 3, 8, 7, 10, 11, 12, 11, 10, 10, 4, 10, 6, 4, 4] ; 
CO2_p_stdev = [22, 20, 20, 19, 18, 19, 17, 14, 13, 14, 14, 16, 18, 18, 17, 16, 17, 14, 14, 13, 11, 11] ; 
O2 = [ 1.9, 5.14, 3.74, 5.05, 14.94, 16.76, 16.22, 24.76, 25.52, 27.25, 31.05, 30.77, 32.38, 28.6, 25.13, 25.86, 26.65, 27.69, 25.12, 21.41, 21.57, 20.72];


% Number of loops %
nsteps = 200 ; 
% Longitude and latitude %
x_lon = 40 ; 
y_lat = 48 ; 

t( 1 ) = 1 ;
time_step = 1 ; %year

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

%Q10 temperature dependence
t25 = 2600 ; 
q_10t = 0.57 ; 

%Intracellular partial pressure 
v = 0.7 ; 

kc25 = 30 ; 
ko25 = 3e4 ; 
q_10c = 2.1 ; 
q_10o = 1.2 ; 
s = 0.015 ; 
alpha = 0.08 ; 
theta = 0.7 ; 
lr_max = 0.75 ; 

%Leaf c:n ratio
CN_leaf = 29 ;

%Tissue respiration rate at 10 degree C
r_tem = 0.055 * 365 ; %0.066 for some PFTs; gC/gN/d -> gC/gN/year
r_bor = 0.066 * 365 ;
r_tro = 0.011 * 365 ; 

%Growth respiraition
R_growth = 0.25 ; 

%Leaf longevity; ranges from 0.5 - 1 depending on type of plant
life_leaf_tem = 0.75 ; 
life_leaf_bor = 0.75 ;
life_leaf_tro = 1 ; 

%Insolation 
ins_present = 150 + 250 .* normpdf( lat', 0, 40 ) ./ normpdf( 0, 0, 40 ) .* ones( x_lon, y_lat ) ;

for a = 1:22

        % Initial time and CO2 setup %
        time = a ; %Ma
        %number between 1-21: 1(-540), 2(-500), 3(-470), 4(-450), 5(-430), 6(-400),
        %7(-370), 8(-340), 9(-300), 10(-280), 11(-260), 12(-245), 13(-220), 14(-200),
        %15(-180), 16(-145), 17(-90), 18(-70), 19(-52), 20(-30), 21(0)

        CO2 = CO2_level(a) ;
        %number between 1-24: 1(10), 2(50), 3(100), 4(160), 5(180), 6(200), 7(240),
        %8(280), 9(360), 10(420), 11(560), 12(700), 13(840), 14(1120), 15(1400), 16(1680), 
        %17(2240), 18(2800), 19(3360), 20(4200), 21(5600), 22(7000), 23(14000),
        %24(28000)

        CO2_ppm = CO2_data( CO2 ) ; 
        %gives the actual CO2 ppm for caculations; ppm = micromol

        %land
        land( : , : ) = double( land_data( : , : , time ) ) ; 
        land( land == 0 ) = NaN ; 

        %temperature for given time and CO2
        tmp_avg( : , : ) = temp_data( : , : , [ CO2 ] , [ time ] ) .* land( : , : ) ;

        %runoff for given time and CO2
        runoff_re = runoff_data( : , : , [ CO2 ] , [ time ] ) .* land( : , : ) ;
        
        % Ice = no land %
        tmp_avg(tmp_avg < -10) = NaN ; 


        %Ambient partial pressure
        pO2 = O2( a ) * 1000 ; %Pa
        pCO2 = CO2_ppm ; %ppmv = micromol

        tf( : , : ) = t25 * ( q_10t .^ ( ( tmp_avg( : , : ) - 25 ) * 0.1 ) ) ;
        tstar( : , : ) = pO2 ./ ( 2 * tf ) ; 
        pi = v * pCO2  ;  
        kc( : , : ) = kc25 * ( q_10c .^ ( ( tmp_avg( : , : ) - 25 ) * 0.1 ) ) ;
        ko( : , : ) = ko25 * ( q_10o .^ ( ( tmp_avg( : , : ) - 25 ) * 0.1 ) ) ;
        c2( : , : ) = ( pi - tstar( : , : ) ) ./ ( pi + kc .* ( 1 + ( pO2 ./ ko ) ) ) ; 

        ftemp_tem( : , : ) = normpdf( tmp_avg( : , : ), 15, 20 ) ;
        ftemp_bor( : , : ) = normpdf( tmp_avg( : , : ), 0, 20 ) ;
        ftemp_tro( : , : ) = normpdf( tmp_avg( : , : ), 27, 7 ) ; 
        c1_bor( : , : ) = alpha * ftemp_bor .* ( ( pi - tstar( : , : ) ) ./ ( pi + 2 * tstar( : , : ) ) ) ; 
        c1_tem( : , : ) = alpha * ftemp_tem .* ( ( pi - tstar( : , : ) ) ./ ( pi + 2 * tstar( : , : ) ) ) ; 
        c1_tro( : , : ) = alpha * ftemp_tro .* ( ( pi - tstar( : , : ) ) ./ ( pi + 2 * tstar( : , : ) ) ) ;

        sigma( : , : ) = ( 1 - ( ( c2( : , : ) - s ) ./ ( c2( : , : ) - theta * s ) ) ) .^ 0.5 ;
        
        %Change in insolation overtime
        ins( : , : ) = ins_present - (ins_present * 4.6/100 * (abs(time_data(a))/570)) ; 
        
        % Maintenance respiration % 

        %Arrhenius equation; temperature function for respiration
        g_T ( : , : ) = exp( 308.56 .* ( ( 1 / 56.02) - ( 1 ./ ( tmp_avg( : , : ) + 46.02 ) ) ) ) ;

        %%% Water presence %%%

        %converting actual runoff into a sigmoidal curve; water presence ranked
        %from 0-1
        test_water = runoff_re ; 
        test_water(test_water == 0 | isnan(test_water)) = NaN ; 
        test_water = circshift( test_water, -6 ) ; 
        water_stress = 1 - (1 ./ ( 1 + exp(0.005 .* (test_water - 450)))) ;

        %%% PHOTOSYNTHESIS %%%
        photosynth_tem( : , : ) = 10 * 365 * ins( : , : ) .* ( c1_tem( : , : ) ./ c2( : , : ) ) .* ( c2( : , : ) - ( 2 * theta - 1 ) * s - 2 .* ( c2( : , : ) - theta * s ) .* sigma( : , : ) ) .* water_stress( : , : ) ; 
        photosynth_bor( : , : ) = 10 * 365 * ins( : , : ) .* ( c1_bor( : , : ) ./ c2( : , : ) ) .* ( c2( : , : ) - ( 2 * theta - 1 ) * s - 2 .* ( c2( : , : ) - theta * s ) .* sigma( : , : ) ) .* water_stress( : , : ) ; 
        photosynth_tro( : , : ) = 10 * 365 * ins( : , : ) .* ( c1_tro( : , : ) ./ c2( : , : ) ) .* ( c2( : , : ) - ( 2 * theta - 1 ) * s - 2 .* ( c2( : , : ) - theta * s ) .* sigma( : , : ) ) .* water_stress( : , : ) ; 
        % gC/ m^2/ year

        %%% CARBON IN LEAF %%%
        C_leaf_tem( : , : , 1 ) = lr_max .* photosynth_tem( : , : ) ;
        C_leaf_bor( : , : , 1 ) = lr_max .* photosynth_bor( : , : ) ;
        C_leaf_tro( : , : , 1 ) = lr_max .* photosynth_tro( : , : ) ;
        % gC/ m^2/ year

        %%% BIOMASS %%%
        biomass_tem( : , : , 1 ) = 2.5e4 .* land( : , : ) ; 
        biomass_tem( biomass_tem == 0 ) = NaN ; 
        biomass_bor( : , : , 1 ) = 2.5e4 .* land( : , : ) ; 
        biomass_bor( biomass_bor == 0 ) = NaN ; 
        biomass_tro( : , : , 1 ) = 2.5e4 .* land( : , : ) ; 
        biomass_tro( biomass_tro == 0 ) = NaN ; 
        % gC/m^2/year

       %Biomass Calculation%
       for n = 1 : nsteps
    
            %%% Leaf respiration %%%
            R_leaf_bor( :, :, n ) = r_bor * ( C_leaf_bor( :, :, n ) / CN_leaf ) .* g_T( : , : ) ;
            R_leaf_bor(R_leaf_bor < 0) = 0 ; 
            R_leaf_tem ( :, :, n ) = r_tem * ( C_leaf_tem( :, :, n ) / CN_leaf ) .* g_T( : , : ) ;
            R_leaf_tem(R_leaf_tem < 0) = 0 ; 
            R_leaf_tro ( :, :, n ) = r_tro * ( C_leaf_tro( :, :, n ) / CN_leaf ) .* g_T( : , : ) ;
            R_leaf_tro(R_leaf_tro < 0) = 0 ; 
        
            %%% NPP %%
            NPP_bor ( :, :, n ) = ( 1 - R_growth ) .* ( photosynth_bor( : , : ) - R_leaf_bor( :, :, n ) ) ;
            NPP_tem ( :, :, n ) = ( 1 - R_growth ) .* ( photosynth_tem( :, : ) - R_leaf_tem( :, :, n ) ) ; 
            NPP_tro ( :, :, n ) = ( 1 - R_growth ) .* ( photosynth_tro( : , : ) - R_leaf_tro( :, :, n ) ) ; 
            %gC/m2/year
        
            %%% Carbon in leaf allocation %%%
            C_leaf_bor( :, :, n+1 ) = ( C_leaf_bor( :, :, n ) .* ( 1 - life_leaf_bor ) ) + ( lr_max .* NPP_bor( :, :, n ) )  ;
            C_leaf_tem( :,:, n+1 ) = ( C_leaf_tem( :, :, n ) .* ( 1 - life_leaf_tem ) ) + ( lr_max .* NPP_tem( :, :, n ) )  ;
            C_leaf_tro( :, :, n+1 ) = ( C_leaf_tro( :, :, n ) .* ( 1 - life_leaf_tro ) ) + ( lr_max .* NPP_tro( :, :, n ) )  ;
        
            %%% Biomass %%%
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
        
        %Saving temp info 
        temp_end( :, :, a ) = tmp_avg ; 
        runoff_end( :, :, a ) = runoff_re ; 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Competition     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Save results     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%collects necessary end result data
% save('FLORA_Paleorun', 'CO2_m_stdev', 'CO2_p_stdev', 'final_biomass', 'biome', 'temp_end', 'runoff_end', 'final_biomass_bor', 'final_biomass_tem', 'final_biomass_tro', 'CO2_level', 'O2', 'time', 'x_lon', 'y_lat', 'final_NPP' )



