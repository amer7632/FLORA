%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Figures & Analysis    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Please download M_PROJ from: 
%https://www.eoas.ubc.ca/~rich/map.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                     FLORA_validation                    %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load FLORA_validation.mat
load FLORA_present_day_data.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COLOURS %%%

biomass_colour = [ 199 234 229
                   128 205 193
                   53 151 143 
                   1 102 94
                   0 60 48 ]./ 234 ; 
compare_colour = [ 84 48 5 
                   140 81 10 
                   191 129 45
                   223 194 125
                   246 232 195
                   245 245 245
                   199 234 229 
                   128 205 193 
                   53 151 143 
                   1 102 94 
                   0 60 48 ] ./ 246 ;

%%% FIGURE 2 %%%
%Biomass comparison

figure
m_proj( 'robinson', 'longitutde' , [ -180 180 ], 'latitude' , [ -90 90 ] )

a = subplot( 3, 1, 1 );
m_pcolor( lon_cru, lat_cru, flipud( log10( c_biomass ) ) )
m_grid
colorbar
colormap( a, biomass_colour ) 
caxis([ 2 4.5 ])
colorbar
title('Actual biomass')

b = subplot( 3, 1, 2 );
m_pcolor( lon_cru, lat_cru, flipud( log10( final_biomass ) ) )
m_grid
colorbar
colormap( b, biomass_colour ) 
caxis([ 2 4.5 ])
colorbar
title('Modelled biomass')

c = subplot( 3, 1, 3 );
m_pcolor( lon_cru, lat_cru, flipud( compare_carbon ) )
m_grid
colorbar
colormap( c, compare_colour ) 
caxis([ -2e4 2e4 ])
colorbar
title('Biomass comparison')

%%% FIGURE 3 %%%
%LAT/LON COMPARISON

figure
subplot(2,1,1)
plot( lon_cru, nansum( flipud( final_biomass ) ) )
hold on
plot( lon_cru, nansum( c_biomass ) )
ylabel( 'g C/ m^{2}' )
title ( 'Longitude' )

subplot(2,1,2)
plot( lat_cru, nansum( flipud( final_biomass ) , 2 ) )
hold on
plot( lat_cru, nansum( flipud( c_biomass ) , 2 ) ) 
ylabel( 'g C/ m^{2}' ) 
title( 'Latitude' ) 
legend( 'Model data', 'Actual data' )
hold off

figure
a1( : , 1 ) =  log10(c_biomass( : )); 
a1( : , 2 ) =  log10(final_biomass( : )) ;
a2 = rmmissing( a1 ) ;
a2( any( isinf( a2 ), 2 ), : ) = [] ; 
coef1 = polyfit( a2( : , 1 ), a2( : , 2 ), 1)
y1 = polyval( coef1, a2( : , 1 ) ) ;
scatter1 = scatter(a2(:,2), a2(:,1),'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerEdgeAlpha', 0.05) ;
hold on
plot( a2( : , 1 ), y1, 'LineWidth', 1.5)

%R-squared for log scale biomass
mdl= fitlm( a2( : , 1 ), a2( : , 2 ) ) ;
%R-squared for raw biomass
mdl2 = fitlm( c_biomass( : ), final_biomass( : ) ) ; 
hold off
xlim([1.5 4.5])
ylim([1.5 4.5])
box on
ylabel('log_{10}(Measured biomass)')
xlabel('log_{10}(Modelled biomass)')

clear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                     FLORA_PALEORUN                      %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load FLORA_Paleorun.mat
load paleo2021_data1.mat
load area.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Colours
biomass_colour = [  223 194 125
                   140 81 10
                   191 127 45
                   246 232 195
                   199 234 229
                   128 205 193
                   53 151 143 
                   1 102 94
                   0 60 48 ] ./ 246 ; 

%%% FIGURE 4 %%%
%Paleobiomass
%Land area was converted to have -1e5 just to illustrate land area with no
%biomass coloured in light brown.

figure
m_proj( 'robinson', 'longitutde' , [ 0 352.5 ] , 'latitude' , [ -87 87 ] )

subplot( 3, 2, 1 )
% 475 Ma
land = double( land_data( :, :, 3 ) ) ; 
land( land == 0 ) = NaN ;
land( land == 1 ) = -1e5 ; 
m_pcolor( lon, lat, circshift( land, 22, 2 ) ) 
hold on
m_pcolor( lon, lat, circshift( log10( final_biomass( :, :, 3) ), 22, 2 ) )
hold off
m_grid( 'fontsize', 0.1 ) 
colormap( biomass_colour )
caxis([2 4.25])
title( 'Ordovician: 470 Ma' )

subplot( 3, 2, 2 )
% 400 Ma
land = double( land_data( :, :, 6 ) ) ; 
land( land == 0 ) = NaN ;
land( land == 1 ) = -1e5 ; 
m_pcolor( lon, lat, land ) 
hold on
m_pcolor( lon, lat, log10( final_biomass( :, :, 6) ) )
hold off
m_grid( 'fontsize', 0.1 ) 
caxis([2 4.25])
title( 'Devonian: 400 Ma' )

subplot( 3, 2, 3 )
% 300 Ma
land = double( land_data( :, :, 9 ) ) ; 
land( land == 0 ) = NaN ;
land( land == 1 ) = -1e5 ; 
m_pcolor( lon, lat, land ) 
hold on
m_pcolor( lon, lat, log10( final_biomass( :, :, 9) ) )
hold off
m_grid( 'fontsize', 0.1 ) 
caxis([2 4.25])
title( 'Carboniferous: 300 Ma' )

subplot( 3, 2, 4 )
% 200 Ma
land = double( land_data( :, :, 14 ) ) ; 
land( land == 0 ) = NaN ;
land( land == 1 ) = -1e5 ; 
m_pcolor( lon, lat, circshift( land, 22, 2 ) )  
hold on
m_pcolor( lon, lat, circshift( log10( final_biomass( :, :, 14) ), 22, 2 ) )
hold off
m_grid( 'fontsize', 0.1 ) 
caxis([2 4.25])
title( 'Jurassic: 200 Ma' )

subplot( 3, 2, 5 )
% 90 Ma
land = double( land_data( :, :, 17 ) ) ; 
land( land == 0 ) = NaN ;
land( land == 1 ) = -1e5 ; 
m_pcolor( lon, lat, circshift( land, 22, 2 ) ) 
hold on
m_pcolor( lon, lat, circshift( log10( final_biomass( :, :, 17) ), 22, 2 ) )
hold off
m_grid( 'fontsize', 0.1 ) 
caxis([2 4.25])
title( 'Cretaceous: 90 Ma' )

subplot( 3, 2, 6 )
% 0 Ma
land = double( land_data( :, :, 22 ) ) ; 
land( land == 0 ) = NaN ;
land( land == 1 ) = -1e5 ; 
m_pcolor( lon, lat, circshift( land, 22, 2 ) ) 
hold on
m_pcolor( lon, lat, circshift( log10( final_biomass( :, :, 22) ), 22, 2 ) )
hold off
m_grid( 'fontsize', 0.1 ) 
caxis([2 4.25])
title( 'Modern: 0 Ma' )
a = colorbar ;
ylabel(a,'log10(Biomass)','FontSize',7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 5 made using data below
%Calculating averages

%latitudinal weights
weights = sum(area,2) ./ sum(sum(area)) ;

for i = 1:22
    %temperature for given time and CO2 (weighted)
    tmp_avg = nanmean(temp_data( : , : , CO2_level(i) ,i) ,2) .* weights  ;
    tmp_m_stdev = nanmean(temp_data(:,:, CO2_m_stdev(i), i),2) .* weights  ; 
    tmp_p_stdev= nanmean(temp_data(:,:,CO2_p_stdev(i), i),2) .* weights  ;
    
    %runoff for given time and CO2 (weighted)
    land = double(land_data(:,:,i)) ;
    land(land == 0) = NaN ; 
    averages(i,1) = nansum(tmp_avg) ;
    averages(i,2) = nansum(tmp_m_stdev) ;
    averages(i,3) = nansum(tmp_p_stdev) ;
    averages(i,4) = sum(nansum(runoff_data(:,:,CO2_level(i),i) .* land )) ; 
    averages(i,5) = sum(nansum(runoff_data(:,:,CO2_m_stdev(i),i) .* land )) ;
    averages(i,6) = sum(nansum(runoff_data(:,:,CO2_p_stdev(i),i) .* land )) ; 
    averages(i,7) = nansum(nansum(area .* land)) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Linear regression analysis

%Calculating r2 values

%Temp vs Biomass


for i = 1:22
    subplot(6,4,i)
    a1 = temp_data(:,:,CO2_level(i),i).* land_data(:,:,i) ; 
    a2 = final_biomass(:,:,i) .* area ; 
    a1= a1(:) ; 
    a2 = a2(:) ; 
    for b = 1:length(a1)
        if isnan(a2(b)) == 1
        a1(b) = -1e5 ; 
        a2(b) = -1e5 ;
        else
        end
    end
    a1(a1 == -1e5 ) = [] ; 
    a2(a2 == -1e5 ) = [] ;

    coef1 = polyfit(a1,a2,1) ;
    y1 = polyval(coef1,a1) ;
    plot(a1, a2, '.') 
    xlim([ -10 35 ])
%     hold on
%     plot(a1, y1)
    mdl= fitlm(a1,a2) ;
    r_sq(i) = mdl.Rsquared ;
end
temp_r = table2array(struct2table(r_sq)) ; 

%Runoff vs Biomass

figure
for i = 1:22
    subplot(6,4,i)
    a1 = runoff_data(:,:,CO2_level(i),i).* land_data(:,:,i) ;
    a2 = final_biomass(:,:,i) .* area ;
    a1= a1(:) ; 
    a2 = a2(:) ; 
    for b = 1:length(a1)
        if isnan(a2(b)) == 1
        a1(b) = -1e5 ; 
        a2(b) = -1e5 ;
        else
        end
    end
    a1(a1 == -1e5 ) = [] ; 
    a2(a2 == -1e5 ) = [] ;
    coef1 = polyfit(a1,a2,1) ; 
    y1 = polyval(coef1,a1) ;
    plot(a1, a2, '.') 
    xlim([ 0 3000 ])
%     hold on
%     plot(a1, y1)
    mdl= fitlm(a1,a2) ;
    r_sq(i) = mdl.Rsquared ;
end
runoff_r = table2array(struct2table(r_sq)) ; 

%Regression with 22 timepoints
for i = 1 : 22
    sum_biomass(i) = sum(nansum(final_biomass(:,:,i) .* area))/1000 ; %kg C
end
coef1 = polyfit(averages(:,1),sum_biomass(:), 1) ; 
y1 = polyval(coef1, averages(:,1)) ; 

figure
plot(averages(:,1), sum_biomass, '.')
xlabel('Temp (^{o}C)')
ylabel('Biomass (kg C)')
hold on
plot(averages(:,1), y1)
hold off
mdl = fitlm(averages(:,1), sum_biomass) 

figure
plot(averages(:,4), sum_biomass, '.')
coef1 = polyfit(averages(:,4), sum_biomass(:),1) ; 
y1 = polyval(coef1, averages(:,4)) ; 
hold on
plot(averages(:,4), y1)
xlabel('Runoff (mm)')
ylabel('Biomass (kg C)')
mdl = fitlm(averages(:,4), sum_biomass)

%NPP over time
for i = 1 : 22
    sum_NPP(i) = sum(nansum(final_NPP(:,:,i) .* area)) ;
end
sum_NPP = sum_NPP/sum_NPP(end) ; 
figure
plot(time_data, sum_NPP)
xlabel('Time (Ma)')
ylabel('Relative NPP')
ylim([0 3])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Biomes         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

area_biomass = final_biomass .* area ; 
for a = 1:time
    for i = 1 : x_lon
        for j = 1 : y_lat
            if biome( i, j, a ) == 1 %Temperate
                btem(i,j,a) = area_biomass(i,j,a) ; 
            elseif biome( i, j, a ) == 2 %Boreal
                bbor(i,j,a) = area_biomass(i,j,a) ; 
            elseif biome(i,j,a) == 3 %Tropical
                btro(i,j,a) = area_biomass(i,j,a) ; 
            else
            end
        end
    end
    
    tem_sum(a) = sum(nansum(btem(:,:,a))) ;
    bor_sum(a) = sum(nansum(bbor(:,:,a))) ; 
    tro_sum(a) = sum(nansum(btro(:,:,a))) ;
    
end



