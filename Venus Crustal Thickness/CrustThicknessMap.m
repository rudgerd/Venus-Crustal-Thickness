%% Path setup
clear all; close all;
% set up path to function folders
current_folder = pwd; func = append(current_folder,'/functions'); 
path(func,path)

%% Venus Crustal thickness (One-Layer model and Uniform Crustal Density)
H = load('VenusData/VenusTopo719.shape');           % Loading Topography data
Clm = load('VenusData/shgj180u.txt');               % Loading Gravity data
Clm = [0 0 0 0 0 0; Clm];
R = .6051000000000000E+07;                          % Reference Radius

% Input Variables
dW = 15e3;                                          % Mean crustal thickness 
rhoc = 2800;                                        % Crustal Density
rhom = 3300;                                        % Mantle Density
lmax = 80;                                          % SH Filter                                      
nmax = 8;
dres = 1;

% Topography Map
H_trunc = H(1:addmup(lmax),1:4);
[tmap,lon,lat] = plm2xyz(H_trunc,dres);

%% Venus Crustal thickness (Two-Layer model and Uniform Crustal Density)
% Input Variables for Two-Layer Model
drhocm = 3000;                                      % Mantle-Core Density Contrast
lw = 70;                                             % SH Filter Crust-Mantle Boundary
lp = 40;                                            % SH Filter Mantle Interfaces
ViscProf = 'isoviscous';                            % Viscosity Profile
dM = 500e3;                                         % Depth of the mantle bottom
planet = 'Venus';                                      

% Two Layer crustal thickness model w/ uniform crustal density
W = TwoLayer(dW,dM,planet,rhoc,rhom,drhocm,lw,lp,lmax,ViscProf); 

[W_map,lon,lat] = plm2xyz(W,dres);%,c11cmn
TwoLayerUniDensity = (tmap-W_map)*10^-3;

figure
imagesc(lon,lat,TwoLayerUniDensity)
set(gca,'YDir','normal')
a = colorbar;
a.Label.String = 'km';
hold on;
contour(lon,lat,TwoLayerUniDensity,4,'LineWidth',1,'LineColor','k');
xlabel('Longitude')
ylabel('Latitude')
title(' Two-Layer Crustal Thickness Model (Uniform Crustal Density)')

%% Make Rho Map
Latitudes = 90:-dres:-90;
Longitudes = 0:dres:360;

% Map out tessera terrain regions
rho = load('VenusData/Tessera_Density_Map.txt');
[rho_map,lon,lat] = plm2xyz(rho,dres);
figure
imagesc(Longitudes,Latitudes,rho_map)
b = colorbar;
b.Label.String = 'kg/m^3';
title('Variable Density Map')
set(gca,'YDir','normal')

%% Two-Layer Crustal Thickness (Lower Crustal Density at Tesserae)

W = TwoLayer(dW,dM,planet,rho,rhom,drhocm,lw,lp,lmax,ViscProf); 
[W_map,lon,lat] = plm2xyz(W,dres);%,c11cmn

TwoLayerVaryDensity = (tmap-W_map)*10^-3;

figure
imagesc(lon,lat,TwoLayerVaryDensity)
set(gca,'YDir','normal')
c = colorbar;
c.Label.String = 'km';
hold on;
contour(lon,lat,TwoLayerVaryDensity,4,'LineWidth',1,'LineColor','k');
xlabel('Longitude')
ylabel('Latitude')
title ('Two-Layer Crustal Thickness (Lower Crustal Density at Tesserae)');

%% The difference between a uniform and varying density modeled

dif = TwoLayerVaryDensity - TwoLayerUniDensity;
figure
imagesc(lon,lat,dif)
set(gca,'YDir','normal')
d = colorbar;
d.Label.String = 'km';
hold on;
contour(lon,lat,dif,4,'LineWidth',1,'LineColor','k');
xlabel('Longitude')
ylabel('Latitude')
title ('Crustal Thickness Difference (Varying Density Model - Uniform Density Model)'); 
