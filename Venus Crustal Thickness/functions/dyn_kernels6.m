function [Z,G,H,C,K,S,vcTH,vrTH]=dyn_kernels6(deg,Bp,Visc,Te,planet,SurfBC,rhoR)
% function [Z,G,H,C,K,S,vcTH,vrTH]=DYN_KERNELS6(deg,Bp,Visc,Te,planet,SurfBC)
%
% Calculates the flow in a viscous, self-gravitating sphere due to a 
% harmonic load, and returns a set of response kernels related to the
% resulting dynamic topography and geoid.  This version is based on a 6x6
% system of equations (see various papers by Hager & O'Connell.
% The effects of bending and membrane stresses are included using
% the formalism of Turcotte et al. 1981.
% Uses the function 'fralmanac.m' to retrieve planetary parameter values
% for Earth and 'ReturnPropagator.m' to calculate propagator matrices.
%
% INPUT:
% 
% deg      Desired spherical harmonic degree (can be a vector of degrees).
% Bp       Radius of the loading interface normalized by the planetary
%          radius (should be less than 1 and greater than the normalized
%          CMB radius).  Can also be a vector of loading radii, in which 
%          case the radii should be sorted in decreasing order.
% Visc     A two/three-columned matrix in which the first column contains
%          normalized interface radii, the second column contains
%          viscosities normalized by the reference viscosity, and the third
%          column is density normalized by reference density. OR: 'PREM',
%          'isoviscous' [default]
%          (Note that none of the returned kernels depend on the absolute 
%          value of the reference viscosity, only relative viscosity
%          contrasts)
% Te       Elastic thickness, in units of km [default: 0]
% planet   String containing the name of the desired planet, from which 
%          physical constants are assumed, e.g. 'Earth', 'Mercury', 'Mars',
%          'Venus', or 'Moon' [default: 'Earth']
%          Alternatively, planet is a vector of parameters:
%          
% SurfBC   0 No-slip surface boundary condition
%          1 Free-slip surface boundary condition [default]
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
% dyn_kernels('demo1') - Reproduces figures 9.21(a,d) of Hager and Clayton
%                        1989
% dyn_kernels('demo2') - Reproduces figures C1(b,c) of Katzman et al. 1998
% dyn_kernels('demo3') - Reproduces figures 6(a,b) of Herrick and Phillips
%                        1992
% 
% Created by Peter James pjames at alum dot mit dot edu, 6/10/12
% Modified 10/27/14 to add the kernel S
% Modified 12/3/14 to add the kernel C

defval('deg',2:20)
defval('Bp',1:-.01:.8)
defval('Visc','isoviscous')
defval('Te',0)
defval('planet','Moon')
defval('rhom',3200) % Mantle density
defval('drho_cm',6000) % Core-mantle density contrast
defval('mu0',10^19)
defval('SurfBC',0)
defval('poisson',.25) % Poisson's ratio, dimensionless
defval('E',10^11) % Young's modulus, units of Pa
% defval('rhoR',3000)
if isempty(strfind(deg(:)','demo'))

if ischar(planet)
    if strcmp(planet,'Venus')
        R = 6051.880596e3;
        RC = 3000e3;
        gr = 8.87; rhom=3400;
    elseif strcmp(planet,'Mercury')
        R = 2440e3;
        RC = 2020e3;
        gr = 3.7;
        rhom=3400;
        drho_cm=4000;
        if min(Bp)<(R-400e3)/R % Assume CMB compensation 
            disp('NOTICE: Compensating at the CMB in dyn_kernels6')
            RC=800e3; % Inner core radius (somewhat arbitrary)
        end
    elseif strcmp(planet,'Mars')
        R = 3396e3;
        RC = 1500e3;
        gr = 3.71;
    elseif strcmp(planet,'Earth')
        R=6371e3;
        RC=3483e3;
        gr=9.8;
    elseif strcmp(planet,'Moon')
        R = 1737e3;
        RC = 200e3;
        gr = 1.622;
        rhom=3200;
    else
        warning('Planet name not recognized')
    end
else
    if length(planet<5)
        error('Input for ''planet'' should be a string or a vector')
    else
        disp('input parameter ''planet'' is a vector:');
        disp('[R RC gr rho_m drho_cm ]');
        R=planet(1);
        RC=planet(2);
        gr=planet(3);
        rhom=planet(4);
        drho_cm=planet(5);
    end
end
% disp(['Radius of ',planet,' (km): ',num2str(R/1000,6)])
% disp(['CMB of ',planet,' (km): ',num2str(C/1000,6)])
% disp(['Gravity on ',planet,' (m/s^2): ',num2str(gr,3)])

% Normalized core radius
RCp=RC/R;

% Convert to meters
Te=Te*1000;

if ischar(Visc)
    rhoR = rhom; % Set the reference density equal to mantle density
    if strcmp(Visc,'PREM')
        Visc = ...
         [(R-670e3)/R, 10, 1;...
         (R-400e3)/R, 1, 1;...
         (R-60e3)/R, 1/30, 1;...
         R/R, 100, 1];
    elseif strcmp(Visc,'isoviscous')
        Visc = [R/R, 1, 1];
    elseif strcmp(Visc,'200km lid')
        Visc = [(R-100e3)/R, 1, 1;...
         R/R, 100, 1];
    elseif strcmp(Visc,'300km lid')
        Visc = [(R-300e3)/R, 1, 1;...
         R/R, 10, 1];
    elseif strcmp(Visc,'500km lid')
        Visc = [(R-500e3)/R, 1;...
         R/R, 10, 1];
    elseif strcmp(Visc,'300km tapered')
        Visc = [(R-300e3)/R, 1, 1;...
            (R-200e3)/R, 10, 1;...
            (R-100e3)/R, 100, 1;...
         R/R, 100, 1];
%     elseif strcmp(Visc,'300km_log10_100km')
%         ct=0;
%         for z=0:10:300
%             ct=ct+1;
%             Visc(ct,1:2) = [(R-(300-z)*10^3)/R, 10^(z/100)];
%         end
%     elseif strcmp(Visc,'500km_log10_100km')
%         ct=0;
%         for z=0:10:500
%             ct=ct+1;
%             Visc(ct,1:2) = [(R-(500-z)*10^3)/R, 10^(z/100)];
%         end
    elseif strcmp(Visc,'Exponential_50km')
        clear Visc
        ct=0;
        folddepth=50;
        maxd=floor((R-RC)/1000/folddepth)*folddepth;
        depths=0:10:maxd;
        for z=depths
            ct=ct+1;
            Visc(ct,1:2) = [(R-(maxd-z)*10^3)/R, 10^(z/folddepth)];
            Visc(ct,3) = 1;
        end
    elseif strcmp(Visc,'Exponential_200km')
        clear Visc
        ct=0;
        folddepth=200;
        maxd=floor((R-RC)/1000/folddepth)*folddepth;
        depths=0:10:maxd;
        for z=depths
            ct=ct+1;
            Visc(ct,1:2) = [(R-(maxd-z)*10^3)/R, 10^(z/folddepth)];
            Visc(ct,3) = 1;
        end
    elseif strcmp(Visc,'Exponential_300km')
        clear Visc
        ct=0;
        folddepth=300;
        maxd=floor((R-RC)/1000/folddepth)*folddepth;
        depths=0:10:maxd;
        for z=depths
            ct=ct+1;
            Visc(ct,1:2) = [(R-(maxd-z)*10^3)/R, 10^(z/folddepth)];
            Visc(ct,3) = 1;
        end
    elseif strcmp(Visc,'Exponential_500km')
        clear Visc
        ct=0;
        folddepth=500;
        maxd=floor((R-RC)/1000/folddepth)*folddepth;
        depths=0:10:maxd;
        for z=depths
            ct=ct+1;
            Visc(ct,1:2) = [(R-(maxd-z)*10^3)/R, 10^(z/folddepth)];
            Visc(ct,3) = 1;
        end
%         maxd
%         disp(max(Visc(:,2)))
%         disp(min(Visc(:,2)))
%         figure
%         plot(flipud(Visc(:,2))/10,depths)
%         ylabel('Depth')
%         xlabel('Normalized Viscosity')
%         box off
%         set(gca,'YDir','reverse','XAxisLocation','top','YTIck', 0:200:1000)
%         ylim([0 1500]
% %         xlim([0 1000])

    elseif strcmp(Visc,'500km tapered')
        Visc = [(R-500e3)/R, 1, 1;...
            (R-400e3)/R, 10, 1;...
            (R-300e3)/R, 100, 1;...
            (R-200e3)/R, 1000, 1;...
            (R-100e3)/R, 10000, 1;...
         R/R, 10000, 1];
    elseif strcmp(Visc,'Lunar')
        clear Visc
        %%
        rho_m=3359;
        rho_c=2800;
        rhoR=3000;
        maxd=800e3;
        dinterval=100e3;
        depths=vpa(0:dinterval:maxd); depths=fliplr(depths);
        Visc=0;
        for i=1:length(depths)
            if depths(i)>100e3
                Visc(i,1:3) = [(R-depths(i))/R, ...
                    1+2*(800e3-depths(i))/700e3, rho_m/rhoR];
            else
                Visc(i,1:3) = [(R-depths(i))/R, ...
                    3+27*(100e3-depths(i))/100e3, rho_m/rhoR];
            end
            if depths(i)<40e3
                Visc(i,3)=rho_c/rhoR;
            end
        end
    else
        warning('String ''Visc'' not understood')
    end
    
end

try
    rho_core = Visc(1,3)*rhoR + drho_cm; % Lower mantle density plus contrast
catch
    rhoR = rhom; % Set the reference density equal to mantle density
    rho_core = Visc(1,3)*rhoR + drho_cm; % Lower mantle density plus contrast
end
rho_cstar = rho_core/rhoR;

Grav=6.6743*10^(-11); % Gravitational constant, units of m^3/kg/s^2
% Assume the gravitational acceleration throughout the interior is 
% approximately equal to the surface acceleration.
gc=gr;
gb=gr;
% Unitary driving load (kg/m^2)
massload=1;

ndegs=length(deg);
ndepths=length(Bp);
% Initialize output matrices
Z=zeros(ndepths,ndegs);
G=zeros(ndepths,ndegs);
H=zeros(ndepths,ndegs);
C=zeros(ndepths,ndegs);
K=zeros(ndepths,ndegs);
S=zeros(ndepths,ndegs); % Ratio of shear potential function over tau_rr
vrTH=zeros(ndepths,ndegs);
vcTH=zeros(ndepths,ndegs);

% L=1:lmax;
for i=1:length(Bp)
RB=Bp(i)*R;
lcount=0;
tic
for ll=1:length(deg)
    l=deg(ll);
    lcount=lcount+1;
    massload_star=massload/(R*rhom);   
    [PBR,PCR]=ReturnPropagator6(Bp(i),RCp,l,Visc);
    % Construct the system of equations
    
    % Nondimensionalization:
    % V{star} = V * mu_0 / (g_r*rho_m*R^2)
    % dr{star} = dr / R
    % Tau{star} = Tau / (g_r*rho_m*R)
    % massload{star} = massload / (rho_m*R)
    gammaR = 4*pi*Grav*R/(gr*(2*l+1));
    gammaC = 4*pi*Grav*RC/(gc*(2*l+1));
    gammaB = 4*pi*Grav*RB/(gb*(2*l+1));

    k1 = (-l^3*(l+1)^3+4*l^2*(l+1)^2)/(-l*(l+1)+1-poisson);
    k2 = (-l*(l+1)+2)/(-l*(l+1)+1-poisson);
    Elastic = E*Te^3/(12*(1-poisson^2)*rhom*gr*(R-Te/2)^4) * k1 ...
        + E*Te/(rhom*gr*(R-Te/2)^2) * k2;
    A=zeros(4,4);
    b=zeros(4,1);
    
    for ii=1:4
        di2=0; di3=0; di4=0;
        if ii==2; di2=1; elseif ii==3; di3=1; elseif ii==4; di4=1; end;
        A(ii,1) = PCR(ii,2);
        A(ii,2) = (PCR(ii,3) + gammaC*rhoR * ...
            (-PCR(ii,3)*rho_cstar + PCR(ii,5) - (l+1)*PCR(ii,6))) * ...
            drho_cm/rhom*gc/gr*RC/R;
        A(ii,3) = di3*(1 + Elastic) + gammaR*RC/R*rhoR * (RC/R)^l * ...
            (-PCR(ii,3)*rho_cstar + PCR(ii,5) + l*PCR(ii,6));
        if SurfBC==0
            A(ii,4) = -di4;
        else
            A(ii,4) = -di2;
        end
        b(ii) = (-PBR(ii,3) + (2*l+1)*gammaB*rhoR*PBR(ii,6) + ...
            gammaB*RC/RB * (RC/RB)^l * rhoR * ...
            (PCR(ii,3)*rho_cstar - PCR(ii,5) - l*PCR(ii,6)))...
            * RB/R*gb/gr * massload_star;
    end

    % Solution vector components: vcTH (lateral velocity at the core), 
    % perturbation at RC, perturbation at R, Tau_theta at the surface

    % Improve the condition number of the system (
    ScaleFactor=NaN(4,1);
    for ii=1:4
        ScaleFactor(ii)=norm(A(:,ii))*10^15;
        A(:,ii)=A(:,ii)/ScaleFactor(ii);
    end
    
    Xvec=A\b;
    
    Xvec=Xvec./ScaleFactor;
    dC_star=Xvec(2);
    dH_star=Xvec(3);

    % Re-dimensionalize:
    dH = dH_star * R;
    dC = dC_star * R;
    % Core velocities and surface velocities / shear stresses, if you  
    % desire them:
    vcTH = Xvec(1) * gr*rhom*R^2 / mu0;
    if SurfBC==0
        tauTH = Xvec(4) * gr*rhom*R;
        S(i,lcount) = tauTH / (dH*gr*rhom); % or Xvec(4)/Xvec(3);
    else
        vrTH = Xvec(4) * gr*rhom*R^2 / mu0;
    end
%     massload = massload_star * R*rho_m;

    UR = 4*pi*Grav/(2*l+1)*(R*rhom*dH + RC*(RC/R)^(l+1)*drho_cm*dC + ...
        RB*(RB/R)^(l+1)*massload);
    UC = 4*pi*Grav/(2*l+1)*(R*(RC/R)^l*rhom*dH + RC*drho_cm*dC + ...
        RB*(RC/RB)^l*massload);
    N=UR/gr;
    dUB = 4*pi*Grav*RB/(2*l+1)*(RB/R)^(l+1)*massload; 
    dB = 4*pi*Grav*RB/(2*l+1)*massload;
    
    Z(i,lcount)=UR/(gr*dH);  
    G(i,lcount)=UR/dB;
    H(i,lcount)=dH*rhom/massload;
    
    C(i,lcount)=dC*drho_cm/massload;
    K(i,lcount)=UR/dUB;
    
%     % Check the math
%     lhs = [0; 0; -rho_m*gr*(dH-UR/gr) - Elastic*rho_m*gr*dH; 0];
%     if SurfBC==0
%         lhs(4)=-tauTH;
%     else
%         lhs(2)=-vrTH*mu0/R;
%     end
%     rhs1 = [0; mu0/R*vcTH; RC/R*drho_cm*gc*(dC-UC/gc); 0];
%     rhs2 = [0; 0; RB/R*gb*massload; 0];
%     
%     rhs = PCR*rhs1 + PBR*rhs2;
%     
%     difer(rhs-lhs);
% %     if abs(rhs-lhs)/norm(lhs)>.001
%         disp('rhs')
%         disp(rhs)
%         disp('lhs')
%         disp(lhs)
% %     end
    
end
toc

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


elseif strcmp(deg,'demo1')
    R=6371; % in km
    CMB=3483; % in km
    depths=0:10:R-CMB;
    Bp = (R*ones(size(depths))-depths)/R;
    [~,G,H,~]=dyn_kernels6([2,4,8],Bp,'isoviscous',0,'Earth',1);
    figure
    plot(depths,-H(:,1))
    hold on
    plot(depths,-H(:,2),':')
    plot(depths,-H(:,3),'--')

    view(90,90)
    hleg1 = legend('L = 2','L = 4','L = 8');
    title('L = 2,4,8 displacement kernels with depth','FontWeight','bold',...
        'FontSize',12)
    xlabel('Loading depth (km)','FontSize',10)
    ylabel('Surface displacement kernels (dimensionless)','FontSize',10)
    ylim([0 1])
    xlim([0 R-CMB])

    figure
    plot(depths,G(:,1))
    hold on
    plot(depths,G(:,2),':')
    plot(depths,G(:,3),'--')
    
    view(90,90)
    hleg2 = legend('L = 2','L = 4','L = 8');
    title('L = 2,4,8 geoid kernels with depth','FontWeight','bold',...
        'FontSize',12)
    xlabel('Loading depth (km)','FontSize',10)
    ylabel('Geoid kernels (dimensionless)','FontSize',10)
    ylim([-.5 .5])
    xlim([0 R-CMB])

elseif strcmp(deg,'demo2')
    depths=0:10:1500;
    R=6371;
    Bp = (R*ones(size(depths))-depths)/R;
    Visc = ...
         [(R-670)/R, 10, 1;...
         (R-400)/R, 1, 1;...
         (R-60)/R, 1/30, 1;...
         R/R, 100, 1];
    [~,G1,H1,~]=dyn_kernels6(26,Bp,Visc);
    [~,G2,H2,~]=dyn_kernels6(26,Bp,Visc,40);
    [~,G3,H3,~]=dyn_kernels6(26,Bp,Visc,130);
    figure
    plot(depths,G1)
    hold on
    plot(depths,G2,'.')
    plot(depths,G3,'--')

    view(90,90); legend;
    hleg1 = legend('Te = 0','Te = 40','Te = 130');
    title('L = 26 geoid kernels with depth','FontWeight','bold',...
        'FontSize',12)
    xlabel('Loading depth (km)','FontSize',10)
    ylabel('Geoid kernel (dimensionless)','FontSize',10)
    ylim([-.3 .3])
    xlim([0 1500])

    figure
    plot(depths,H1)
    hold on
    plot(depths,H2,'.')
    plot(depths,H3,'--')
    
    view(90,90);
    hleg2 = legend('Te = 0','Te = 40','Te = 130');
    title('L = 26 displacement kernels with depth',...
        'FontWeight','bold','FontSize',12)
    xlabel('Loading depth (km)','FontSize',10)
    ylabel('Surface displacement kernel (dimensionless)','FontSize',10)
    ylim([-1 1])
    xlim([0 1500])
    
elseif strcmp(deg,'demo3')
    deg=2:18;
    R=6051;
    Bp=(R-200)/R;
    [Z1,~,~,K1]=dyn_kernels6(deg,Bp,'isoviscous',0,'Venus');
    Visc2 = [(R-800)/R, 1, 1;...
             R/R, .1, 1];
    [Z2,~,~,K2]=dyn_kernels6(deg,Bp,Visc2,0,'Venus');
    Visc3 = [(R-800)/R, 1, 1;...
             R/R, .01, 1];
    [Z3,~,~,K3]=dyn_kernels6(deg,Bp,Visc3,0,'Venus');
    Visc4 = [(R-200)/R, 1, 1;...
             R/R, .1, 1];
    [Z4,~,~,K4]=dyn_kernels6(deg,Bp,Visc4,0,'Venus');
    Visc5 = [(R-200)/R, 1, 1;...
             R/R, .01, 1];
    [Z5,~,~,K5]=dyn_kernels6(deg,Bp,Visc5,0,'Venus');

    figure
    plot(deg,K1)
    hold on
    plot(deg,K2,'linestyle','--','Color','r')
    plot(deg,K3,'linestyle','--','Color','g')
    plot(deg,K4,'linestyle','--','Color','c')
    plot(deg,K5,'linestyle','--','Color','m')
    hleg1 = legend('isoviscous','Model B','Model C','Model D','Model E');
    title('Potential kernels at 200 km depth',...
        'FontWeight','bold','FontSize',12)
    xlabel('Spherical harmonic degree','FontSize',10)
    ylabel('K','FontSize',10)
    ylim([-1 1])
    xlim([2 18])

    figure
    plot(deg,Z1)
    hold on
    plot(deg,Z2,'linestyle','--','Color','r')
    plot(deg,Z3,'linestyle','--','Color','g')
    plot(deg,Z4,'linestyle','--','Color','c')
    plot(deg,Z5,'linestyle','--','Color','m')
    hleg2 = legend('isoviscous','Model B','Model C','Model D','Model E');
    title('Admittance kernels at 200 km depth',...
        'FontWeight','bold','FontSize',12)
    xlabel('Spherical harmonic degree','FontSize',10)
    ylabel('Z','FontSize',10)
    ylim([-.1 .1])
    xlim([2 18])
% elseif strcmp(deg,'demo4')
% % Try to get fix divergence at high spherical harmonic degrees
% l=61;


end

end

