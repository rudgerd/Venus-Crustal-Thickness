function [Z,G,H,C,K,S,vcTH,vrTH]=...
    dyn_kernels_density(deg,Rvec,Visc,Te,planet,SurfBC,rhoR)
% function [Z,G,H,C,K,S,vcTH,vrTH]=...
% DYN_KERNELS_DENSITY(deg,Bp,Visc,Te,planet,SurfBC)
%
% Calculates dynamic flow kernels for a density anomaly distributed over a
% range of mantle depths.  This function calls dyn_kernels
%
% INPUT:
% 
% deg      Desired spherical harmonic degree (can be a vector of degrees).
% Rvec     A vector of length 2 or 3, where the first element is the
%          shallower end, the second is the deeper end, and the third is
%          the depth increment [NORMALIZED BY R]
% Visc     A two-columned matrix in which the first column contains
%          normalized interface radii and the second column contains
%          viscosities normalized by the reference viscosity. OR: 'PREM',
%          'isoviscous' [default]
%          (Note that none of the returned kernels depend on the absolute 
%          value of the reference viscosity, only relative viscosity
%          contrasts)
% Te       Elastic thickness, in units of km [default: 0]
% planet   String containing the name of the desired planet, from which 
%          physical constants are assumed, e.g. 'Earth', 'Mercury', 'Mars',
%          'Venus', or 'Moon' [default: 'Earth']
% SurfBC   0 No-slip surface boundary condition [default]
%          1 Free-slip surface boundary condition
% 
% OUTPUT:
%
% Z        A non-dimensional (ndepths)x(ndegs) matrix giving the
%          degree-dependent ratios of geoid to topography for each inputted
%          loading depth (the "admittance kernel").
% G        A non-dimensional (ndepths)x(ndegs) matrix giving the "geoid
%          kernel".  This represents the ratio of gravitational potential 
%          for a given mass load, relative to the static potential anomaly 
%          produced by a similar mass anomaly at the surface.
% H        A non-dimensional (ndepths)x(ndegs) matrix giving the normalized
%          surface displacements.
% C        A non-dimensional (ndepths)x(ndegs) matrix giving the normalized
%          CMB displacements.
% K        A non-dimensional (ndepths)x(ndegs) "potential kernel" used by
%          some authors (related to the Geoid kernel G by an
%          upward-continuation factor.
%
% Created by Peter James pjames at alum dot mit dot edu, 5/10/16

% Mostly, use the defaults in dyn_kernels
defval('deg',2:30)
% defval('Rvec',[100e3, 300e3, 1e3])
defval('Visc',[])
defval('Te',[])
defval('planet','Moon')
defval('rho_m',3200)
defval('drho_cm',6000)
defval('mu0',[])
defval('SurfBC',[])
defval('poisson',[])
defval('E',[])
defval('rhoR',3000)
Grav=6.6743*10^(-11); % Gravitational constant, units of m^3/kg/s^2

if ~ischar(deg)

if ischar(planet)
    if strcmp(planet,'Venus')
        R = 6051.880596*10^3;
        RC = 3000*10^3;
        gr = 8.87; 
        rho_m=3400;
    elseif strcmp(planet,'Mercury')
        R = 2440*10^3;
        RC = 2020*10^3;
        gr = 3.7;
        rho_m=3400;
        drho_cm=4000;
    elseif strcmp(planet,'Mars')
        R = 3396*10^3;
        RC = 1500*10^3;
        gr = 3.71;
    elseif strcmp(planet,'Earth')
        R=fralmanac('Radius',planet);
        RC=fralmanac('CMB',planet);
        gr=fralmanac('GravAcc',planet);
    elseif strcmp(planet,'Moon')
        R = 1737*10^3;
        RC = 200*10^3;
        gr = 1.622;
    else
        warning('Planet name not recognized')
    end
else
    warning('Please enter a string for the input parameter ''planet''');
end

if length(Rvec)==2
    Bp=Rvec(2):1000/R:Rvec(1);
else
    Bp=Rvec(2):Rvec(3):Rvec(1);
end

% [Z,G,H,C,K,S,vcTH,vrTH]=dyn_kernels(deg,Bp,Visc,Te,planet,SurfBC);
% Outputs are arrays [ndepths]x[lmax]
[O1,O2,O3,O4,O5,O6,O7,O8]=dyn_kernels6(deg,Bp,Visc,Te,planet,SurfBC,rhoR);

    
Z=zeros(length(Bp),length(deg));
G=Z; K=Z; UR=Z;

massload=1;
RB=zeros(1,length(Bp));

for i=1:length(Bp)
    for il=1:length(deg)
        l=deg(il);
        RB(i)=R*Bp(i);
        dUB = 4*pi*Grav*RB(i)/(2*l+1)*(RB(i)/R)^(l+1)*massload; 
        dB = 4*pi*Grav*RB(i)/(2*l+1)*massload;
        dH = O3(i,il)*massload/rho_m;
    %     dC(i,il) = C*load/drho_cm;
        UR = O2(i,il)*dB;

        Z(i,il)=UR/(gr*dH);  
        G(i,il)=UR/dB;
    %     H(i,il)=dH*rho_m/load;
    %     C(i,il)=dC*drho_cm/load;
        K(i,il)=UR/dUB;    
    end
    problem_index=find(O3(i,2:end)>0,1);
    if ~isempty(problem_index)
        % Uh-oh, the propagator matrix problem is poorly conditioned
        % Extrapolate the exponential decay: H = -A*exp(-b*l)
        disp(['Extrapolating l=' int2str(problem_index+1)])
        j1=floor(problem_index/2);
        if j1<3
            warning('The kernels are irredeemable!')
        end
        logslope=log(abs(O3(i,j1+1)))-log(abs(O3(i,j1)));
        for j=j1+1:size(O3,2)
            O3(i,j) = O3(i,j1+1)*exp((j-j1)*logslope);
            G(i,j) = O3(i,j)*Z(i,j);
        end
    end
end


% Z = mean(Z,1);
% G = mean(G,1);
% H = mean(O3,1);
% C = mean(O4,1);
% K = mean(K,1);
% S = mean(O6,1);
% vcTH = mean(O7,1);
% vrTH = mean(O8,1);

% Weight by depth:
Z = RB.^2 * Z / sum(RB.^2);
G = RB.^2 * G / sum(RB.^2);
H = RB.^2 * O3 / sum(RB.^2);
C = RB.^2 * O4 / sum(RB.^2);
K = RB.^2 * K / sum(RB.^2);
S = RB.^2 * O6 / sum(RB.^2);
vcTH = RB.^2 * O7 / sum(RB.^2);
vrTH = RB.^2 * O8 / sum(RB.^2);

% % Find the integrated geoid and topography from a unitary density
% % anomaly
% figure
% plot(Z,'k')
% hold on
% plot(H,'b')
% plot(G,'g')
% title('Admittance/Displacement/Geoid')
% xlim([2 50]); ylim([-1 1])



elseif strcmp(deg, 'demo1')
    deg = 2:60;
    R = 1737e3;
    Rvec = [R-100e3,R-300e3,1e3]/R;
    Bp1 = (R-100e3)/R;
    Bp2 = (R-200e3)/R;
    Bp3 = (R-300e3)/R;

    Visc='isoviscous';
    Te=0;
    planet='Moon';
    SurfBC=1;
    
    tic
    [Z,G,H,C,K,S,vcTH,vrTH]=...
        dyn_kernels_density(deg,Rvec,Visc,Te,planet,SurfBC);
    toc

    tic
    [Z1,G1,H1,C1,K1,S1,vcTH1,vrTH1]=...
        dyn_kernels(deg,Bp1,Visc,Te,planet,SurfBC);
    toc

    tic
    [Z2,G2,H2,C2,K2,S2,vcTH2,vrTH2]=...
        dyn_kernels(deg,Bp2,Visc,Te,planet,SurfBC);
    toc
    
    tic
    [Z3,G3,H3,C3,K3,S3,vcTH3,vrTH3]=...
        dyn_kernels(deg,Bp3,Visc,Te,planet,SurfBC);
    toc
    
    figure
    plot(Z,'k')
    hold on
    plot(Z1,'b')
    plot(Z2,'g')
    plot(Z3,'r')
    title('Admittance')

    figure
    plot(G,'k')
    hold on
    plot(G1,'b')
    plot(G2,'g')
    plot(G3,'r')
    title('Geoid Kernel')

    figure
    plot(H,'k')
    hold on
    plot(H1,'b')
    plot(H2,'g')
    plot(H3,'r')
    title('Topography Kernel')
    
elseif strcmp(deg,'demo2')
    R=1737.4e3; 
    lmax=60;
    degrees=1:lmax;
    rho=load('DATA/grain_density_310.sh'); 
    rho(:,3:4)=rho(:,3:4)*.88;
    rhoc0=rho(1,3); % Degree zero crustal density
    %     rho=2560;
    rhom0=3220;
%     drhocm=5000;
    planet='Moon';
    Te=0;
%     [~,GM,R] = read_TAB('GRAV/JGGRAIL_900C11A_SHA.TAB');
    GM = 4.902800222140801e+12;
    R=1738000;
    g=GM/R^2; gw=g; gr=g; gc=g;
    dW=34e3; dV=1400e3;
    Rw=R-dW;
    Rm=R-dV;
    depthinterval=10000;
%     Rc = 200*10^3;
    rhoR=3000; % Reference density
    SurfBC=0;
%         rhoc=2560; drho=500;
% ViscProf = [(R-800e3)/R, 1, rhom0/rhoR;...
%     (R-700e3)/R, 1.25, rhom0/rhoR;...
%     (R-600e3)/R, 1.50, rhom0/rhoR;...
%     (R-500e3)/R, 1.75, rhom0/rhoR;...
%     (R-400e3)/R, 2.00, rhom0/rhoR;...
%     (R-300e3)/R, 2.33, rhom0/rhoR;...
%     (R-200e3)/R, 2.66, rhom0/rhoR;...
%     (R-100e3)/R, 3, rhom0/rhoR;...
%     Rw/R, 27*(100e3-(R-Rw))/100e3+3, rhom0/rhoR;...
%     R/R, 30, rhoc0/rhoR];  
% ViscProf='Lunar';
    maxd=800e3;
    dinterval=100e3;
    depths=vpa(0:dinterval:maxd); depths=fliplr(depths);
    ViscProf=vpa(0);
    for i=1:length(depths)
        if depths(i)>100e3
            ViscProf(i,1:3) = vpa([(R-depths(i))/R, ...
                1+2*(800e3-depths(i))/700e3, rhom0/rhoR]);
        else
            ViscProf(i,1:3) = vpa([(R-depths(i))/R, ...
                3+27*(100e3-depths(i))/100e3, rhom0/rhoR]);
        end
        if depths(i)<40e3
            ViscProf(i,3)=vpa(rhoc0/rhoR);
        end
    end
%     ViscProf='Lunar';

    Rvec = [Rw/R, Rm/R, depthinterval/R];
    tic
    [Z,G,H,C] = dyn_kernels_density(degrees,Rvec,ViscProf,Te,planet,SurfBC);
    toc
%     Dmantle_R = drH * gw/gr;
%     Dmantle_Rc = drC * gc/gr;
    
    Bp1 = (R-dW)/R;
    Bp2 = (R-mean([dW dV]))/R;
    Bp3 = (R-dV)/R;

    tic
    [Z1,G1,H1,C1,K1,S1,vcTH1,vrTH1]=...
        dyn_kernels6(degrees,Bp1,ViscProf,Te,planet,SurfBC);
    toc
    
    tic
    [Z2,G2,H2,C2,K2,S2,vcTH2,vrTH2]=...
        dyn_kernels6(degrees,Bp2,ViscProf,Te,planet,SurfBC);
    toc
    
    tic
    [Z3,G3,H3,C3,K3,S3,vcTH3,vrTH3]=...
        dyn_kernels6(degrees,Bp3,ViscProf,Te,planet,SurfBC);
    toc
    
    figure
    plot(Z,'k')
    hold on
    plot(Z1,'b')
    plot(Z2,'g')
    plot(Z3,'r')
    title('Admittance')
    xlim([2 50]); ylim([0 .6])

    figure
    plot(G,'k')
    hold on
    plot(G1,'b')
    plot(G2,'g')
    plot(G3,'r')
    title('Geoid Kernel')
    xlim([2 50]); ylim([-1 1])

    figure
    plot(H,'k')
    hold on
    plot(H1,'b')
    plot(H2,'g')
    plot(H3,'r')
    title('Topography Kernel')
    xlim([2 50]); ylim([-1 1])

else
    warning('ON','Input not recognized')
end

end

