function [W,rhom,Hmantle,Umantle,H,RGH] = ...
    TwoLayer(dW,dM,planet,rhoc,rhom0,drhocm,lw,lp,lmax,ViscProf)
%
% Perform an inversion for a two-layered model of crust & mantle structure. 
% This implements "Model B", i.e., the two interior mass anomalies are 
% associated with (1) relief on the crust-mantle interface and (2) lateral
% variations in mantle density.
% We endeavor to find the mass anomalies that minimize Bouguer gravity 
% misfit while simultaneously minimizing uncompensated topography.
%
% INPUTS: 
% dW - Mean crustal thickness
% dM - Bottom of the upper mantle, in which density is allowed to vary. 
% planet - planet name
% rhoc - crustal density (either a constant, a map, or coefficients)
% rhom0 - mantle density (constant)
% drhocm - density contrast between the mantle and outer core
% lw - critical degree for filtering W
% lp - critical degree for filtering V
% Te_Cl - value of elastic thickness for degree of compensation
% f - amplitude of top loading relative to V
%
% OUTPUTS:
% W - Crust-mantle interface relief
% rhom - Global density Map
% W1 - One-layered solution
% RGHmodel - reproduced value for rhoc*g*H
% Nmodel - reproduced value for the Geoid
% H - Topography, for your convenience
% RGH - Weight of surface topography, for your convenience
% psiF - Mass column anomaly not supported by Moho relief, Mantle mass
%        anomaly, or self-gravitation [kg/m^2]
%
% Last modified 10/10/16 by pjames at alum.mit.edu
% Last modified 04/05/2020 by rhdame@gmail.com
% Input Variables for One-Layer Model

if ~ischar(dW)

defval('rhom0',3200)
defval('drhocm',5000)
defval('dW',40e3)
defval('dV',300e3)
defval('lw',80);
defval('lp',30);
defval('wreduction',.5);
defval('vreduction',.5);
defval('resfactor',2); % Factor by which the xyz2plm grid res increases
defval('planet','Venus')
defval('Te',0)
defval('f',1)
defval('finiteamplitude',0)
defval('depthinterval',10000);
defval('SurfBC',0)
defval('mindegV',3)
defval('MasconMask',0)

if ischar(planet) % planet gives the name of the planet 
    if strcmp(planet,'Mercury')
%         H = load('SHAPE/fitswp164_120.out');
%         H(:,3:4)=H(:,3:4)*10^3;
%         [~,GM,R] = read_Geodyn('GRAV/grvfld.HGM003c');
%         Clm=load('GRAV/HgM_EGU15.txt');
        H = read_TAB('SHAPE/gtmes_150v05_sha.tab');
        H(:,3:4)=H(:,3:4)*10^3;
        [Clm,GM,R] = read_TAB('GRAV/ggmes_100v07_sha.tab');
        U = plm2pot(Clm,R,GM,R,1,'spherical');
        Rc = 2020*10^3;
        defval('ViscProf','200km lid');
        defval('lmax',90);
    elseif strcmp(planet,'Venus')
        H = load('VenusData/VenusTopo719.shape');
        RH = 6051.880596000000*10^3;
        Clm = load('VenusData/shgj180u.txt');
        R = 0.6051000000000000e07;
        GM = 0.3248585920790000e15;
        Clm = [0 0 0 0 0 0; Clm];
        U = plm2pot(Clm,R,GM,R,1,'spherical');
        Rc = 3000*10^3;
        defval('ViscProf','200km lid');
        defval('lmax',90);
    elseif strcmp(planet,'Earth')
        H = load('SHAPE/srtmp360.ret_ellipsoid');
        R = 6368637.46025739;
        GM = 3.986005*10^14;
        Clm = load('GRAV/EGM2008_ZeroTide360'); 
        U = plm2pot(Clm,R,GM,R,1,'spherical');
        Rc=fralmanac('CMB',planet);
        defval('ViscProf','200km lid');
        defval('lmax',90);
    elseif strcmp(planet,'Moon')
        H = load('SHAPE/LOLA1439v2.sh');
        [Clm,GM,R] = read_TAB('GRAV/JGGRAIL_900C11A_SHA.TAB');
        U = plm2pot(Clm,R,GM,R,1,'spherical');
        defval('ViscProf','Lunar');
        Rc = 200*10^3;
        defval('lmax',90);
        MasconMask=1;
    elseif strcmp(planet,'Mars')
        H = load('SHAPE/MarsTopo719.shape');
        H(4,3)=.05*H(4,3);
        H(6,3:4)=.25*H(6,3:4);
        [Clm,GM,R]=read_TAB('GRAV/jgmro_110c_sha.tab');
        Clm(4,3)=.05*Clm(4,3);
        Clm(6,3:4)=.25*Clm(6,3:4);
        U = plm2pot(Clm,R,GM,R,1,'spherical');
        defval('ViscProf','200km lid');
        defval('lmax',90);
        Rc = 1500*10^3;
%         rhoc=2700; drho=500;
    else
        warning('Planetary body name not recognized')
    end
end
% defval('gr',GM/R^2)
% defval('gw',GM/R^2)
% defval('gv',GM/R^2)
% defval('gc',GM/R^2)

% While "R" is the radius at which the gravity field is calculated, "RH" is
% the mean radius of the planetary shape.
if H(1,3)>R/2 % H(1,3) equals RH
    RH = H(1,3);
else % Assume this is a correction to R
    RH = R+H(1,3);
end
% RH=R;
Rw=RH-dW;
Rm=RH-dM;
try
    gr = GravityAcceleration(planet,R);
    gw = GravityAcceleration(planet,Rw);
    gv = GravityAcceleration(planet,Rm);
    gc = GravityAcceleration(planet,Rc);
catch
    disp('Using surface acceleration everywhere')
    gr=GM/R^2;
    gw=gr; gv=gr; gc=gr;
end


Grav=6.67300*10^-11;
defval('mindeg',1)
% dres=180/lmax; %sqrt(lmax*(lmax+1)); % Nyquist
% dres2=dres/resfactor; % Resolution of xyz2plm inputs
dres2=1; % This will be sufficient for all expansions up to lmax=180

degrees=1:lmax; degrees=degrees';

% Truncate at lmax
H=H(1:addmup(lmax),1:4);
U=U(1:addmup(lmax),1:4);
% Remove the lowest degrees
H(1:addmup(mindeg-1),3:4)=zeros(addmup(mindeg-1),2);
U(1:addmup(mindeg-1),3:4)=zeros(addmup(mindeg-1),2);
% Create lat-lon map
Hmap=plm2xyz(H,dres2);
% Umap=plm2xyz(U,dres2);

% Deal with densities, define density maps
% Crustal density
if size(rhoc,2)>1 && size(rhoc,2)<7 % rhoc is lmcosi matrix
    rhoc0=rhoc(1,3); % Degree zero crustal density
    try rhoc=rhoc(1:addmup(lmax),:); % rhoc has more than lmax degrees
    catch % rhoc has fewer than lmax degrees, do nothing
    end
    rhoc_map=plm2xyz(rhoc,dres2);
elseif size(rhoc,2)>6 % rhoc is a map
    rhoc_map=rhoc;
    rhoc=xyz2plm(rhoc,lmax);
    rhoc0=rhoc(1,3);
else  % rhoc is a scalar
    rhoc_map=rhoc*ones(size(Hmap));
    size(rhoc_map)
    rhoc0=rhoc;
    % NOTE: no rhoc_lm defined
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Make an array of viscous properties, if not already inputed

rhoR=3000; % Reference density
if ischar(ViscProf)
    if strcmp(ViscProf,'isoviscous')
        ViscProf = [Rm/RH, 1, rhom0/rhoR;...
            Rw/RH, 1, rhom0/rhoR;...
            RH/RH, 1, rhoc0/rhoR];       
    elseif strcmp(ViscProf,'PREM')
        ViscProf = [(RH-670e3)/RH, 10, rhom0/rhoR;...
            (RH-400e3)/RH, 1, rhom0/rhoR;...
            (RH-60e3)/RH, 1/30, rhom0/rhoR;...
            RH/RH, 100, 1];
    elseif strcmp(ViscProf,'200km lid')
        ViscProf = [Rm/RH, 1, rhom0/rhoR;...
            (RH-200e3)/RH, 1, rhom0/rhoR;...
            Rw/RH, 1, rhom0/rhoR;...
            RH/RH, 1, rhoc0/rhoR];       
    elseif strcmp(ViscProf,'300km lid')
        ViscProf = [Rm/RH, 1, rhom0/rhoR;...
            (RH-300e3)/RH, 1, rhom0/rhoR;...
            Rw/RH, 1, rhom0/rhoR;...
            RH/RH, 1, rhoc0/rhoR];       
    elseif strcmp(ViscProf,'500km lid')
        ViscProf = [Rm/RH, 1, rhom0/rhoR;...
            (RH-500e3)/RH, 1, rhom0/rhoR;...
            Rw/RH, 1, rhom0/rhoR;...
            RH/RH, 1, rhoc0/rhoR];       
    elseif strcmp(ViscProf,'300km tapered')
        ViscProf = [Rm/RH, 1, rhom0/rhoR;...
            (RH-300e3)/RH, 1, rhom0/rhoR;...
            (RH-200e3)/RH, 10, rhom0/rhoR;...
            (RH-100e3)/RH, 100, rhom0/rhoR;...
            Rw/RH, 1000, rhom0/rhoR;...
            RH/RH, 1000, rhoc0/rhoR];  
    elseif strcmp(ViscProf,'Lunar')
        ViscProf = [(RH-800e3)/RH, 1, rhom0/rhoR;...
            (RH-700e3)/RH, 1.286, rhom0/rhoR;...
            (RH-600e3)/RH, 1.571, rhom0/rhoR;...
            (RH-500e3)/RH, 1.857, rhom0/rhoR;...
            (RH-400e3)/RH, 2.143, rhom0/rhoR;...
            (RH-300e3)/RH, 2.429, rhom0/rhoR;...
            (RH-200e3)/RH, 2.714, rhom0/rhoR;...
            (RH-100e3)/RH, 3, rhom0/rhoR;...
            Rw/RH, 27*(100e3-(RH-Rw))/100e3+3, rhom0/rhoR;...
            RH/RH, 30, rhoc0/rhoR];  
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Calculate displacements resulting from all mass anomalies

% First, find filters for the solutions so we can pick lmax2
omega_w=zeros(lmax,1); omega_v=zeros(lmax,1); Gwl=zeros(lmax,1); Gml=Gwl;
% Max degree for 2x2 inversion (less than the final solution)
lmax2=min([lmax 60]); 
for l=1:lmax
    Gwl(l)=4*pi*Grav/gr/(2*l+1) * Rw*(Rw/R)^(l+1);
    rmantle=RH-200e3; %Rw; %100e3;
    Gml(l)=4*pi*Grav/gr/(2*l+1) * rmantle*(rmantle/R)^(l+1);
end
for l=1:lmax
    lam_w = Gwl(lw)^2*(1/wreduction-1);
    omega_w(l) = 1 / (1 + lam_w/Gwl(l)^2);
%     lam_v = Gwl(lp)^2*(1/vreduction-1);
%     omega_v(l) = 1 / (1 + lam_v/Gwl(l)^2);
    lam_v = Gml(lp)^2*(1/vreduction-1);
    omega_v(l) = 1 / (1 + lam_v/Gml(l)^2);
    if lmax2==min([lmax 60]) && omega_v(l)<1e-3
        lmax2=l;
    end
end
% Round lmax2 to the nearest 10
lmax2=round(lmax2/10)*10;
disp('lmax2:')
disp(lmax2)

degrees2=1:lmax2;

% Surface & core displacements from the Moho
option=1;
if option==0 % dynamic flow kernels
    [~,~,drH,drC] = dyn_kernels6(degrees2,Rw/RH,ViscProf,Te,planet,SurfBC,rhoR);
    Dw_R = drH * gw/gr;
    Dw_Rc = drC * gc/gr;
elseif option==1 % Equal pressures (approximately Cartesian equal-mass)
    Dw_R=ones(lmax2,1);
    % Equal pressures from Hemingway & Matsuyama (2017)
    for l=1:lmax2
        Dw_R(l) = -gw/gr;        
    end
    Dw_Rc=zeros(lmax2,1);
elseif option==2
    Dw_R=ones(lmax2,1);
    for l=1:lmax2
        % Two options according to Phillips and Sleep (1985):
        % Eliminate shear stresses: PS85 = (Rw/R)^2
        % Minimize deviatoric stress: PS85 = 1-dW/R*3*l*(l+1)/(2*l*(l+1)-1);
        PS85 = 1 - dW/RH * 3*l*(l+1)/(2*l*(l+1)-1);
        Dw_R(l) = -PS85 * gw/gr;
    end
    Dw_Rc=zeros(lmax2,1);
end

% Surface & core displacements from crustal density
if option==0 % dynamic flow kernels
    Rvec = [RH/RH, Rw/RH, depthinterval/RH];
    [~,~,drH,drC] = dyn_kernels_density(degrees2,Rvec,ViscProf,Te,planet,SurfBC,rhoR);
    Dcrust_R = drH * gw/gr;
    Dcrust_Rc = drC * gc/gr;
elseif option==1 % Equal pressures (approximately Cartesian equal-mass)
    Dcrust_R=ones(lmax2,1);
    % Equal pressures from Hemingway & Matsuyama (2017)
    for l=1:lmax2
        Dcrust_R(l) = -gw/gr;        
    end
    Dcrust_Rc=zeros(lmax2,1);
end

% Surface & core displacements from mantle density
Rvec = [Rw/RH, Rm/RH, depthinterval/RH];
[~,~,drH,drC] = dyn_kernels_density(degrees2,Rvec,ViscProf,Te,planet,SurfBC,rhoR);
Dmantle_R = drH * gw/gr;
Dmantle_Rc = drC * gc/gr;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Generate masks, if desired

% Optional: create a data window in order to focus on a particular
% feature (in this case, SPA)
SPAwindow=0;
if SPAwindow==1
    % Make a taper
    TH=45;
    L=4;
    dresT=2;
    [G1,V1,EL,EM]=glmalpha(TH,L);
    [G,V]=SortSlep(G1,V1);

    Taper = zeros(size(V));
    for i=1:length(V)
        lmcosi=addmpos(G(:,i),EL,EM);
        map = plm2xyz(lmcosi,dresT);
        lmcosi2 = xyz2plm(map.^2,L);
        if i==1
            Taper=zeros(size(map));
            Tlmcosi=lmcosi2; Tlmcosi(:,3:4)=0;
        end
        Tlmcosi(:,3:4)=Tlmcosi(:,3:4)+lmcosi2(:,3:4)*V(i);
        Taper=Taper+(map.^2)*V(i);
    end
    Tlmcosi(:,3:4) = Tlmcosi(:,3:4)/max(max(Taper));
    Taper = Taper / max(max(Taper));

    SPAlat=-53;
    SPAlon=191;
    SPAcolat=90-SPAlat;
    Tlmcosi=plm2rot(Tlmcosi,0,-SPAcolat,360-SPAlon);
    SPA_taper = plm2xyz(Tlmcosi,dres2);
end

% Optional: Create a mask to remove mascons
if MasconMask==1
    spreadfactor = 4;
    filestr=strcat('DATA/MasconMask_dres1_',int2str(spreadfactor),'x.txt');
    try
        Mascon_taper=load(filestr);
    catch
        Mascons = [...
            37	341.5	11.3
            25.4	18.8	9.2
            16.8	58.4	8.2
           -15.6	35.1	7.3
            -2.5	86.9	7.2
            -20.1	265.2	7.2
            -23.8	320.8	5.9
            -4.6	52	5.9
            25.1	190.6	5.5
            51.2	237.5	5.4
            -49.8	265.4	5.4
            18.35	175.2	5.2
            57.26	82	5.1
            -16	293	4.4
            -36.1	208.3	4.4
            2	231	4.2
            26	148	3.7
            13.4	201.8	3.6
            -5	291.3	3.6
            -55.7	314.8	3.5
            39.5	218	3.5
            4.8	23.4	3.4
            -4.4	202.2	3.3
            -57.3	163.1	3.1
            -81	123	2.8
            58.70	213.9	2.8
            34.2	263	2.6
            5.70	140.90	2.6
            -74.9	133.5	2.5
            -32.8	163.8	2.5 
            -15.8	69.6	2.3
            -57.4	135.1	2.1
            -31.25	112.8	2.1
            -35.4	194	2.0
            -21.20  128.9 1.5
            -47.1   176.2   2.1
            -38.2   179.2   2.3
        ];

        MasconNames = {...
            'Imbrium' ...
            'Serenitatis' ...
            'Crisium' ...
            'Nectaris' ...
            'Smythii' ...
            'Orientale' ...
            'Humorum' ...
            'Fecunditatis' ...
            'Fitzgerald-Jackson' ...
            'Coulomb-Sarton' ...
            'Mendel-Rydberg' ...
            'Freundlich-Sharonov' ...
            'Humboldtianum' ...
            'Cruger-Sirsalis' ...
            'Apollo' ...
            'Hertzsprung' ...     
            'Moscoviense' ...     
            'Dirichlet-Jackson' ...
            'Grimaldi' ...
            'Schiller-Zucchius' ...
            'Fowler-Charlier' ...
            'Lamont' ...     
            'Korolev' ...
            'Poincare' ...
            'Amundsen-Ganswindt' ... 
            'Birkhoff' ...
            'Lorentz' ...
            'Mendeleev' ...
            'Schrodinger' ...
            'Ingenii' ...
            'Balmer-Kapteyn' ...
            'Planck' ...
            'Milne' ...
            'Oppenheimer' ...
            'Tsiolkovskiy' ...
            'VonKarman'...
            'Leibniz'...
            };
        Mlat=Mascons(:,1);
        Mlon=Mascons(:,2);
        Mrad=Mascons(:,3);
        Mrad(6)=Mascons(6,3)*1.5;
        Mrad(16)=Mascons(16,3)*1.5;
        dres2=1;
        Mascon_taper = ones(181,361);%size(Hmap));
        Tlat = 90:-dres2:-90;
        Tlon = 0:dres2:360;
        for i=1:size(Mascons,1)
            for j=1:length(Tlat)
                for k=1:length(Tlon)
                    arclen = distance(Tlat(j),Tlon(k),Mlat(i),Mlon(i));
                    if arclen < Mrad(i)*spreadfactor
                        Mascon_taper(j,k) = Mascon_taper(j,k) * ...
                        (1 - cos(pi/2*arclen/(Mrad(i)*spreadfactor))^2);
                    end
                end
            end
            disp(i)
        end
        save(filestr,'Mascon_taper','-ascii')
       
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Calculate gravity-from-topography and gravity-from-density
npower=6;
if finiteamplitude==0
    npower=1;
end
wat=1; % Potential
if length(rhoc)>1 % Crustal density variations are given
    [Utopo,~] = Topo2Grav_PJ(H,rhoc_map,R,RH,GM,R,wat,npower,lmax,lmax*resfactor);
    Ucrust=Density2Grav(rhoc,R,R,Rw,wat,gr);          % ? Ucrust=Density2Grav(rhoc,R,RH,Rw,wat,gr);
else % Crustal density is constant
    [Utopo,~] = Topo2Grav_PJ(H,rhoc,R,RH,GM,R,wat,npower,lmax,lmax*resfactor);
    [~,~,~,Ucrust]=addmon(lmax);
end

% STEP 3: Calculate topographic stress
try % works if rhoc is spatially varying
    RGH_map = gr*rhoc_map.*Hmap;
    RGH = xyz2plm(RGH_map,lmax);
catch
    RGH = H;
    RGH(:,3:4) = rhoc*gr*H(:,3:4);
end

% STEP 4: Calculate the buoyancy produced by crustal density variations
[~,~,~,RGT]=addmon(lmax);
try % works if rhoc is spatially varying
    for i=2:addmup(lmax)
        l=H(i,1);
        RGT(i,3:4)=gr*(RH-Rw)*(rhoc(i,3:4))*Dcrust_R(l);
    end
catch
    % Do nothing
end

% STEP 5: Perform iterations to find W and rhom in two-layered model
% First, iterate to find the stresses consistent with self-gravitation, 
% then calculate gravity correction associated with finite amplitude relief.

imax = 100;
imax2= 100;

[~,~,~,Ures]=addmon(lmax2);
[~,~,~,RGHres]=addmon(lmax2);
[~,~,~,Umoho_finite]=addmon(lmax2);
[~,~,~,Ucmb_finite]=addmon(lmax2);
[~,~,~,psiW_old]=addmon(lmax2);
[~,~,~,psiM_old]=addmon(lmax2);
[~,~,~,psiC_old]=addmon(lmax2);
[~,~,~,psiW]=addmon(lmax2);
[~,~,~,psiM]=addmon(lmax2);
[~,~,~,psiC]=addmon(lmax2);

dampfactor=1; % Average new solution with this amount of the old solution
for iteration=1:imax
    if iteration==imax
        error('iteration maxed out')
    end    
    for i=2:addmup(lmax2)
        l=H(i,1);
        if length(rhoc)>1
            Ures(i,3:4) = U(i,3:4) - Utopo(i,3:4) - Ucrust(i,3:4) -...
                Umoho_finite(i,3:4) - Ucmb_finite(i,3:4) -...
                4*pi*Grav*R/(2*l+1)*(Rc/R)^(l+2)*Dcrust_Rc(l)*...
                (RH-Rw)*rhoc(i,3:4);
        else
            Ures(i,3:4) = U(i,3:4) - Utopo(i,3:4) - Ucrust(i,3:4) -...
                Umoho_finite(i,3:4) - Ucmb_finite(i,3:4);
        end                        
        RGHres(i,3:4) = RGH(i,3:4) - RGT(i,3:4);

        for j=3:4
            bvec = [RGHres(i,j)/gr; ...
                (2*l+1)/(4*pi*Grav*R)*Ures(i,j)];
            Bmat = zeros(2,2);
            Bmat(1,1) = Dw_R(l);
            Bmat(1,2) = Dmantle_R(l);
            Bmat(2,1) = (Rw/R)^(l+2) + (Rc/R)^(l+2)*Dw_Rc(l);
            Bmat(2,2) = (R/(Rw-Rm)/(l+3))*...
                ((Rw/R)^(l+3)-(Rm/R)^(l+3)) + ...
                (Rc/R)^(l+2)*Dmantle_Rc(l);
            X = Bmat\bvec;
            psiW(i,j) = X(1);
            psiM(i,j) = X(2);
        end
    end
    if mindegV>0
        psiM(1:addmup(mindegV-1),3:4)=0;
    end
    % Optional: apply a window to psiV
    if MasconMask==1
        l_lowpass = 0;
        psiM_HP = psiM; 
        psiM_HP(1:addmup(l_lowpass),3:4)=0;
        psiM_LP = psiM; 
        psiM_LP(addmup(l_lowpass)+1:end,3:4)=0;
        if l_lowpass>0
            psiMmap_LP=plm2xyz(psiM_LP,dres2);
        end
        psiMmap_HP=plm2xyz(psiM_HP,dres2);
        if SPAwindow == 1
            psiM=xyz2plm(psiMmap_LP.*SPA_taper + ...
                psiMmap_HP.*SPA_taper.*Mascon_taper,lmax2);
        else
            if l_lowpass>0
                psiM=xyz2plm(psiMmap_LP + ...
                    psiMmap_HP.*Mascon_taper,lmax2);                
            else
                psiM = xyz2plm(psiMmap_HP.*Mascon_taper,lmax2);
            end
        end
    elseif SPAwindow==1
        psiM=xyz2plm(psiVmap.*SPA_taper,lmax2);
    end

    % After masking, apply filter to the solutions
    for i=2:addmup(lmax2)
        l=H(i,1);
        psiM(i,3:4) = omega_v(l) * psiM(i,3:4);
        psiW(i,3:4) = omega_v(l) * psiW(i,3:4);
    end

    for i=2:addmup(lmax2)
        psiC(i,3:4) = Dw_Rc(l)*psiW(i,3:4) + Dmantle_Rc(l)*psiM(i,3:4);
        if length(rhoc)>1
            psiC(i,3:4) = psiC(i,3:4) + Dcrust_Rc(l)*(RH-Rw)*rhoc(i,3:4);
        end
    end
    
    % Dampen the iterated solutions
    psiW(:,3:4)=(psiW(:,3:4)+dampfactor*psiW_old(:,3:4))/(dampfactor+1);
    psiM(:,3:4)=(psiM(:,3:4)+dampfactor*psiM_old(:,3:4))/(dampfactor+1);
    psiC(:,3:4)=(psiC(:,3:4)+dampfactor*psiC_old(:,3:4))/(dampfactor+1);
    [dicW,CF_W]=CheckConvergence(psiW,psiW_old,1e-6);
    [dicM,CF_M]=CheckConvergence(psiM,psiM_old,1e-6);
    if dicW==2 || dicM==2
       error('ON','Solution diverging.  Consider reducing lp and lw')
    end
    if dicW==1 && dicM==1
        disp(['Convergence at iteration = ' num2str(iteration)])
        disp(['Convergence factor for W = ' num2str(CF_W)])
        break
    end
    psiW_old=psiW; psiM_old=psiM; psiC_old=psiC;
	% Finite amplitude corrections
    if iteration<imax
        npower=4;
        wat=1; % Potential anomaly
        % Find mantle density
%         rhom=psiM; rhom(:,3:4)=psiM(:,3:4)/(Rw-Rm); rhom(1,3)=rhom0;
        % NOTE: The interior is composed of cones, not columns
        % Volume of cone with unit area at the surface:
        Vcone=RH/3*((Rw/RH)^3-(Rm/RH)^3);
        rhom=psiM; rhom(:,3:4)=psiM(:,3:4)/Vcone; rhom(1,3)=rhom0;
        % Calculate the new drho
        rhom_map=plm2xyz(rhom,dres2);
%         figure; imagesc(rhom_map); colorbar
        drho_map = rhom_map-rhoc_map;
        if min(min(drho_map))<0
            warning('ON','mantle density is less than crustal density');
        end
        % We solved for psiW, but we need to extract just the relief W
        psiW_map=plm2xyz(psiW,dres2);
        Wmap=psiW_map./drho_map;
        W=xyz2plm(Wmap,lmax2); W(1,3)=0;
        [~,Umoho_finite] = Topo2Grav_PJ(W,drho_map,R,Rw,GM,R,...
            wat,npower,lmax2,lmax2*resfactor);
        C=psiC; C(:,3:4)=psiC(:,3:4)/drhocm;
        [Ucmb,Ucmb_finite] = Topo2Grav_PJ(C,drhocm,R,Rc,GM,R,...
            wat,npower,lmax2,lmax2*resfactor);
    end
end
[~,~,~,Umantle]=addmon(lmax2);
% [~,~,~,Hmantle]=addmon(lmax2);
% Missing the displacement kernel
Hmantle_map = Vcone*(rhom_map-rhom0)./rhoc_map; 
Hmantle=xyz2plm(Hmantle_map,lmax2);
for i=2:addmup(lmax2)
    l=Umantle(i,1);
    Umantle(i,3:4) = 4*pi*Grav*R^2/(2*l+1)/(l+3) * ...
        ((Rw/R)^(l+3) - (Rm/R)^(l+3)) * rhom(i,3:4);
    Hmantle(i,3:4) = Hmantle(i,3:4)*Dmantle_R(l);    
end

        % % Plot power
        % psiMpower=calculate_power(psiM);
        % psiWpower=calculate_power(psiW);
        % figure
        % semilogy(psiMpower,'r')
        % hold on
        % psiMpower2=psiMpower;
        % for i=1:length(psiMpower); psiMpower2(i)=psiMpower(i)/omega_v(i)^2; end
        % semilogy(psiMpower2,'r')
        % semilogy(psiWpower)
        % % ylim([min(psiMpower2) max(psiMpower)])
        % disp(lmax2)

        %figure; plot(omega_v)

% STEP 6: Re-solve for the crust-mantle relief by removing contributions
% of the filtered mantle interface and CMB from the geoid.

[~,~,~,psiW]=addmon(lmax);  % Re-initialize with maximum degree lmax
[~,~,~,Ures]=addmon(lmax);
[~,~,~,Umoho_finite]=addmon(lmax);
for i=2:addmup(lmax)
    l=H(i,1);
        Ures(i,3:4) = U(i,3:4) - Utopo(i,3:4) - Ucrust(i,3:4);
    if l<=lmax2
        Ures(i,3:4) = Ures(i,3:4) - Umantle(i,3:4) - Ucmb(i,3:4);
    end
end
RGHmodel=Umantle;
% Iterate to find finite amplitude crustal thickness
for iteration2=1:imax2
    iteration2
    if iteration2==imax2
        error('iteration2 maxed out')
    end
    for i=2:addmup(lmax)
        l=psiW(i,1);
        psiW(i,3:4) = omega_w(l) * (Ures(i,3:4) - Umoho_finite(i,3:4)) * ...
            (R/Rw)^(l+2)*(2*l+1)/(4*pi*Grav*R);
    end
    % Dampen the iterated solutions
    if iteration2>1
       psiW(:,3:4)=(psiW(:,3:4)+dampfactor*psiW_old(:,3:4))/(dampfactor+1);
        [dicW,CF_W]=CheckConvergence(psiW,psiW_old,1e-6);
        if dicW==1
            disp(['Convergence at iteration2 = ' num2str(iteration2)])
            disp(['Convergence factor for W = ' num2str(CF_W)])
            break
        elseif dicW==2
            error('ON','Solution diverging.  Consider reducing lw')
            break
        end
    end
    psiW_old=psiW;
    % We solved for psiW, but we need to extract just the relief W
    psiW_map=plm2xyz(psiW,dres2);
    Wmap=psiW_map./drho_map;
    W=xyz2plm(Wmap,lmax); W(1,3)=0;
    if iteration2<imax2
        [~,Umoho_finite] = Topo2Grav_PJ(W,drho_map,R,Rw,GM,R,...
            wat,npower,lmax,lmax*resfactor);
    end
end
% Add the radii back in:
W(1,3)
W(1,3)=Rw;
rhom(1,3)=rhom0;

% % STEP 7: Perform calculations for the output parameters
% 
% if ~exist('W','var')
%     psiW_map=plm2xyz(psiW,dres2);
%     Wmap=psiW_map./drho_map;
%     W=xyz2plm(Wmap,lmax);
% end
% 
% if ~exist('V','var')
%     psiV_map=plm2xyz(psiV,dres2);
%     Vmap=psiV_map./drhom_map;
%     V=xyz2plm(Vmap,lmax);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(dW,'demo1')
    lmax=120;
    rho=load('DATA/grain_density_310.sh'); rho(:,3:4)=rho(:,3:4)*.88;
%     rho=2560;
    rhom0=3200;
    drhocm=5000;
    Visc = 'isoviscous';
[Wspa,rhom] = ...
    TwoLayerB(43e3,600e3,'Moon',rho,rhom0,drhocm,80,15,lmax,Visc);

    plotmap(Wspa,[],1); title('W'); colorbar
    plotmap(rhom,[],1); title('rhom'); colorbar
    
    SPAlat=-53;
    SPAlon=191;
    SPAcolat=90-SPAlat;
    rhom_rot=plm2rot(rhom,SPAlon-360,SPAcolat-90,0);
    plotmap(rhom_rot,[],1); title('mantle density rotated'); colorbar
    Wspa1=plm2rot(Wspa,SPAlon-360,SPAcolat-90,0);
    plotmap(Wspa1,[],1); title('Wrotated'); colorbar
    
    %%
    lmax=90;
    Mascon_taper=load('DATA/MasconMask_dres1_4x.txt');
    Mascon = xyz2plm(Mascon_taper,lmax);
%     SPAlat=-53;
%     SPAlon=191;
%     SPAcolat=90-SPAlat;
%     Mascon_rot=plm2rot(Mascon,SPAlon-360,SPAcolat-90,0);

    plotmap(Mascon,[],1); colorbar;
%%
    disp('Done')
    
%     lat1=94; lon1=214;
%     lat2=150; lon2=220;
%     dres=180/lmax;
%     LAT=0:dres:180;
%     LON=0:dres:360;
%     Vmap=plm2xyz(rhom,1);
%     Vmap(LAT==lat2,LON==lon2)
%     Vmap(LAT==lat1,LON==lon1)
%     max(max(Vmap))
%     min(min(Vmap))
%     
%     IsostaticMisfit=RGH;
%     IsostaticMisfit(:,3:4)=RGH(:,3:4)-RGHmodel(:,3:4);
%     
%     plotmap(IsostaticMisfit); title('Isostatic Misfit'); colorbar
%     SMM=calculate_power(IsostaticMisfit);
%     SRR=calculate_power(RGH);
% 
%     figure
%     semilogy(SMM,'r'); title('Isostatic Misfit Power')
%     hold on
%     semilogy(SRR,'k');
% %     CTspa=Hspa;
% %     CTspa(:,3:4)=Hspa(:,3:4)-Wspa(:,3:4);
% % 
% %     CTmap=plm2xyz(CTspa,1);
% %     figure; imagesc(CTmap); colorbar
% 
% %     CTspa1=Hspa;
% %     CTspa1(:,3:4)=Hspa(:,3:4)-W1spa(:,3:4);
% % 
% %     CTmap1=plm2xyz(CTspa1,1);
% %     figure; imagesc(CTmap1); colorbar
% % 
% %     Vmap=plm2xyz(Vspa,1);
% %     figure; imagesc(Vmap); colorbar
% 
%     Wiso=Hspa; Wiso(:,3:4)=Wiso(:,3:4)*-2560/500;
%     
%     
%     
%     Hmap=plm2xyz(Hspa,1);
%     Wmap=plm2xyz(Wspa,1);
%     Wisomap=plm2xyz(Wiso,1);
%     
%     Wmap=Wmap-40e3*ones(size(Wmap));
%     Wisomap=Wisomap-40e3*ones(size(Wisomap));
% 
%     Spp=calculate_power(Vspa);
%     Sww=calculate_power(Wspa);
%     figure
%     semilogy(Spp)
%     hold on
%     semilogy(Sww)
%     
%     Latitudes=[-42 -54 -68];
%     Longitudes=[192];
% %     for i=1:length(Longitudes)
% %         figure
% %         hold on
% %         try
% %             plot((-57:57)*30.3,[Hmap(90:181,Longitudes(i)); ...
% %                 flipud(Hmap(158:180,Longitudes(i)-180))]/1000,'k')
% %             plot((-57:57)*30.3,[Wmap(90:181,Longitudes(i)); ...
% %                 flipud(Wmap(158:180,Longitudes(i)-180))]/1000,'r')
% %             plot((-57:57)*30.3,[Wbfmap(90:181,Longitudes(i)); ...
% %                 flipud(Wbfmap(158:180,Longitudes(i)-180))]/1000,'b')
% %             plot((-57:57)*30.3,[W1map(90:181,Longitudes(i)); ...
% %                 flipud(W1map(158:180,Longitudes(i)-180))]/1000,'k')
% %             plot((-57:57)*30.3,[Wisomap(90:181,Longitudes(i)); ...
% %                 flipud(Wisomap(158:180,Longitudes(i)-180))]/1000,':k')
% %             title(['Longitude ' int2str(Longitudes(i))])
% %         catch
% %             plot((-57:57)*30.3,[Hmap(90:181,Longitudes(i)); ...
% %                 flipud(Hmap(158:180,Longitudes(i)+180))]/1000,'k')
% %             plot((-57:57)*30.3,[Wmap(90:181,Longitudes(i)); ...
% %                 flipud(Wmap(158:180,Longitudes(i)+180))]/1000,'r')
% %             plot((-57:57)*30.3,[Wbfmap(90:181,Longitudes(i)); ...
% %                 flipud(Wbfmap(158:180,Longitudes(i)+180))]/1000,'b')
% %             plot((-57:57)*30.3,[W1map(90:181,Longitudes(i)); ...
% %                 flipud(W1map(158:180,Longitudes(i)+180))]/1000,'k')
% %             plot((-57:57)*30.3,[Wisomap(90:181,Longitudes(i)); ...
% %                 flipud(Wisomap(158:180,Longitudes(i)+180))]/1000,':k')
% %             title(['Longitude ' int2str(Longitudes(i))])
% %         end
% %     end
% %     for j=1:length(Latitudes)
% %         figure
% %         hold on
% %         plot((-110:130)*30.3*cos(Latitudes(j)*pi/180),...
% %             Hmap(90-Latitudes(j),81:321)/1000,'k')
% %         plot((-110:130)*30.3*cos(Latitudes(j)*pi/180),...
% %             Wmap(90-Latitudes(j),81:321)/1000,'r')
% %         plot((-110:130)*30.3*cos(Latitudes(j)*pi/180),...
% %             Wbfmap(90-Latitudes(j),81:321)/1000,'b')
% %         plot((-110:130)*30.3*cos(Latitudes(j)*pi/180),...
% %             W1map(90-Latitudes(j),81:321)/1000,'k')
% %         plot((-110:130)*30.3*cos(Latitudes(j)*pi/180),...
% %             Wisomap(90-Latitudes(j),81:321)/1000,':k')
% %         title(['Latitude ' int2str(Latitudes(j))])
% %     end
    
else
    warning('Input W not recognized')
end

end

