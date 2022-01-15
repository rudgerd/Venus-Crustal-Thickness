function [y1,y2,y3,y4,y5,y6]=dyn_flow6(deg,rp,Rvp,Visc,Te,planet,SurfBC)
% function [y1,y2,y3,y4,y5,y6]=DYN_FLOW6(deg,rp,Rvp,Visc,Te,planet,SurfBC)
%
% Calculates the flow in a viscous, self-gravitating sphere due to a 
% harmonic load, and returns components of the dynamic flow state vector
% This version is based on a 6x6 system of equations (see various papers by 
% Hager & O'Connell).
% The effects of bending and membrane stresses are included using
% the formalism of Turcotte et al. 1981.
% Uses the function 'fralmanac.m' to retrieve planetary parameter values
% for Earth and 'ReturnPropagator6.m' to calculate propagator matrices.
%
% INPUT:
% 
% deg      Desired spherical harmonic degree.
% rp       Receiver radii at which the state vector is measured, normalized
%          by the planetary radius ("p" stands for "prime")
% Rv       Radius of the loading interface normalized by the planetary
%          radius (should be less than 1 and greater than the normalized
%          CMB radius).
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
%          'Venus', or 'Moon' [default: 'Moon']
%          Alternatively, planet is a vector of parameters:
%          
% SurfBC   0 No-slip surface boundary condition
%          1 Free-slip surface boundary condition [default]
% 
% OUTPUT:
%
% u1       Normalized radial velocity (vr_lm)
% u2       Normalized poloidal velocity (phi_lm)
% u3       Normalized normal stress (sigma_rr_lm)
% u4       Normalized shear stress (S_lm)
% u5       Normalized gravitational potential (U_lm)
% u6       Normalized radial derivative of gravitational potential (gr_lm)
% 
%
% dyn_kernels('demo1') - Reproduces figures 9.21(a,d) of Hager and Clayton
%                        1989
% dyn_kernels('demo2') - Reproduces figures C1(b,c) of Katzman et al. 1998
% dyn_kernels('demo3') - Reproduces figures 6(a,b) of Herrick and Phillips
%                        1992
% 
% Created by Peter James pjames at alum dot mit dot edu, 11/19/17

defval('deg',2:20)
defval('Rv',1:-.01:.8)
defval('Visc','isoviscous')
defval('Te',0)
defval('planet','Moon')
defval('rho_m',3200) % Mantle density
defval('rho_m',2800) % Crust density
defval('drho_cm',3000) % Core-mantle density contrast
defval('SurfBC',1)
defval('poisson',.25) % Poisson's ratio, dimensionless
defval('E',10^11) % Young's modulus, units of Pa
defval('precise',0)
muR = 10^19;

if isempty(strfind(deg(:)','demo'))

if ischar(planet)
    if strcmp(planet,'Venus')
        R = 6051.880596e3;
        RC = 3000e3;
        gr = 8.87; rho_m=3400;
        drho_cm=4000;
    elseif strcmp(planet,'Mercury')
        R = 2440e3;
        RC = 2020e3;
        gr = 3.7;
        rho_m=3400;
        drho_cm=4000;
    elseif strcmp(planet,'Mars')
        R = 3396e3;
        RC = 1500e3;
        gr = 3.71;
        rho_m=3300;
        drho_cm=4000;
    elseif strcmp(planet,'Earth')
        R=6371e3;
        RC=3483e3;
        gr=9.8;
        drho_cm=4000;
    elseif strcmp(planet,'Moon')
        R = 1737.2e3;
        RC = 373e3;
        gr = 1.622;
        rho_m=3359;
        rho_c=2800;
        drho_cm=5871-rho_m;
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
        rho_m=planet(4);
        drho_cm=planet(5);
    end
end

if ischar(Visc)
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
    elseif strcmp(Visc,'500km tapered')
        Visc = [(R-500e3)/R, 1, 1;...
            (R-400e3)/R, 10, 1;...
            (R-300e3)/R, 100, 1;...
            (R-200e3)/R, 1000, 1;...
            (R-100e3)/R, 10000, 1;...
         R/R, 10000, 1];
    elseif strcmp(Visc,'LunarVisc')
        % For a reference viscosity of 10^19:
        % Equals 10^27 Pa s at the surface, decreases by a factor of ten 
        % every 100 km down to a minimum viscosity of 10^19 below 800 km
        % depth
        clear Visc
        folddepth=100e3;
        maxd=800e3;
        dinterval=10e3;
        depths=vpa(0:dinterval:maxd); depths=fliplr(depths);
        for i=1:length(depths)
            Visc(i,1:2) = [(R-depths(i))/R, 10^((maxd-depths(i))/folddepth)];
            Visc(i,1:2) = [(R-depths(i))/R, 1];
            if depths(i)<40e3
                Visc(i,3) = 1;%2734/rho_m;
            else
                Visc(i,3) = 1;
            end
        end
        
%         Visc = [0.5, 1, 1;...
%         0.7, 1, 1;...
%         1, 1, 1];
    elseif strcmp(Visc,'Lunar')
        clear Visc
        maxd=800e3;
        dinterval=10e3;
        depths=vpa(0:dinterval:maxd); depths=fliplr(depths);
        rhoR=rho_m;
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
%         Visc = [(R-800e3)/R, 1, rhom0/rhoR;...
%             (R-700e3)/R, 1.25, rhom0/rhoR;...
%             (R-600e3)/R, 1.50, rhom0/rhoR;...
%             (R-500e3)/R, 1.75, rhom0/rhoR;...
%             (R-400e3)/R, 2.00, rhom0/rhoR;...
%             (R-300e3)/R, 2.33, rhom0/rhoR;...
%             (R-200e3)/R, 2.66, rhom0/rhoR;...
%             (R-100e3)/R, 3, rhom0/rhoR;...
%             Rw/R, 27*(100e3-(R-Rw))/100e3+3, rhom0/rhoR;...
%             R/R, 30, rhoc0/rhoR]; 
    else
        warning('String ''Visc'' not understood')
    end
    
end

% Normalized core radius
RCp=RC/R;

rho_R = rho_m; % Set the reference density equal to the mantle density
rho_core = Visc(1,3)*rho_R + drho_cm; % Lower mantle density plus contrast
rho_cstar = rho_core/rho_R;

Grav=6.6743*10^(-11); % Gravitational constant, units of m^3/kg/s^2
% Assume the gravitational acceleration throughout the interior is 
% approximately equal to the surface acceleration.
gc=gr;
gb=gr;
% Unitary driving load (kg/m^2)

precise=0;
if precise==1
    load=vpa(1);
else
    load=1;
end

ndegs=length(deg);
ndepths=length(rp);
nloads=length(Rvp);
% Initialize output matrices
y1=zeros(ndepths,ndegs,nloads);
y2=y1; y3=y1; y4=y1; y5=y1; y6=y1; 

for k=1:nloads
    Rv=Rvp(k)*R;
for ll=1:length(deg)
    if precise==1
        l=vpa(deg(ll));
    else
        l=deg(ll);
    end
    load_star=load/(R*rho_m);   
    
    [PBR,PCR]=ReturnPropagator6(Rvp(k),RCp,l,Visc);
    
    % Construct the system of equations
    
    % Nondimensionalization:
    % V{star} = V * mu_0 / (g_r*rho_m*R^2)
    % dr{star} = dr / R
    % Tau{star} = Tau / (g_r*rho_m*R)
    % load{star} = load / (rho_m*R)
    gammaR = 4*pi*Grav*R/(gr*(2*l+1));
    gammaC = 4*pi*Grav*RC/(gc*(2*l+1));
    gammaB = 4*pi*Grav*Rv/(gb*(2*l+1));

    k1 = (-l^3*(l+1)^3+4*l^2*(l+1)^2)/(-l*(l+1)+1-poisson);
    k2 = (-l*(l+1)+2)/(-l*(l+1)+1-poisson);
    Elastic = E*Te^3/(12*(1-poisson^2)*rho_m*gr*(R-Te/2)^4) * k1 ...
        + E*Te/(rho_m*gr*(R-Te/2)^2) * k2;
    if precise==1
        A=vpa(zeros(4,4));
        b=vpa(zeros(4,1));
    else
        A=zeros(4,4);
        b=zeros(4,1);
    end
    
    for ii=1:4
        di2=0; di3=0; di4=0;
        if ii==2; di2=1; elseif ii==3; di3=1; elseif ii==4; di4=1; end
        A(ii,1) = PCR(ii,2);
        A(ii,2) = (PCR(ii,3) + gammaC*rho_R * ...
            (-PCR(ii,3)*rho_cstar + PCR(ii,5) - (l+1)*PCR(ii,6))) * ...
            drho_cm/rho_m*gc/gr*RC/R;
        A(ii,3) = di3*(1 + Elastic) + gammaR*RC/R*rho_R * (RC/R)^l * ...
            (-PCR(ii,3)*rho_cstar + PCR(ii,5) + l*PCR(ii,6));
        if SurfBC==0
            A(ii,4) = -di4;
        else
            A(ii,4) = -di2;
        end
        b(ii) = (-PBR(ii,3) + (2*l+1)*gammaB*rho_R*PBR(ii,6) + ...
            gammaB*RC/Rv * (RC/Rv)^l * rho_R * ...
            (PCR(ii,3)*rho_cstar - PCR(ii,5) - l*PCR(ii,6)))...
            * Rv/R*gb/gr * load_star;
    end

    % Solution vector components: vcTH (lateral velocity at the core), 
    % perturbation at RC, perturbation at R, Tau_theta at the surface

%     % Improve the condition number of the system (
%     ScaleFactor=NaN(4,1);
%     for ii=1:4
%         ScaleFactor(ii)=1;%norm(A(:,ii))*10^0;
%         A(:,ii)=A(:,ii)/ScaleFactor(ii);
%     end
    
    Xvec=vpa(A)\vpa(b);
    
%     Xvec=Xvec./ScaleFactor;
    vcTH_star=Xvec(1);
    dC_star=Xvec(2);
    dH_star=Xvec(3);
    if SurfBC==0
        tauTH_star=Xvec(4);
        vrTH_star=0;
    elseif SurfBC==1
        tauTH_star=0;
        vrTH_star=Xvec(4);
    end        
    % Re-dimensionalize:
    dH = dH_star * R;
    dC = dC_star * R;
    % Core velocities and surface velocities / shear stresses, if you  
    % desire them:
    vcTH = vcTH_star * gr*rho_m*R^2 / muR;
    vrTH = vrTH_star * gr*rho_m*R^2 / muR;
    tauTH = tauTH_star * gr*rho_m*R;
    
    U_R = 4*pi*Grav/(2*l+1)*(R*rho_m*dH + RC*(RC/R)^(l+1)*drho_cm*dC ...
        + Rv*(Rv/R)^(l+1)*load);
    U_C = 4*pi*Grav/(2*l+1)*(R*(RC/R)^l*rho_m*dH + RC*drho_cm*dC ...
        + Rv*(RC/Rv)^l*load);
    dUdr_C = 4*pi*Grav/(2*l+1)*(l*(RC/R)^(l-1)*rho_m*dH ...
        - (l+1)*drho_cm*dC ...
        + l*(RC/Rv)^(l-1)*load);
    dUdr_R = 4*pi*Grav/(2*l+1)*(l*rho_m*dH ...
        - (l+1)*drho_cm*dC*(RC/R)^(l+2) ...
        - (l+1)*(Rv/R)^(l+2)*load);
    
    % Reconstruct the state vector at the CMB
    u_C(1,1) = 0;
    u_C(2,1) = vcTH;
    u_C(3,1) = RC/muR * (dC*drho_cm*gc - rho_core*U_C);
    u_C(4,1) = 0;
    u_C(5,1) = RC*rho_R/muR * U_C;
    u_C(6,1) = RC^2*rho_R/muR * dUdr_C;
        
    du_B(1,1) = 0;
    du_B(2,1) = 0;
    du_B(3,1) = Rv/muR*gb * load;
    du_B(4,1) = 0;
    du_B(5,1) = 0;
    du_B(6,1) = -4*pi*Grav*Rv^2*rho_R/muR * load;
    
    u_R(1,1) = 0;
    u_R(3,1) = -R/muR*rho_R*gr*dH;
    u_R(5,1) = R*rho_R/muR * U_R;
    u_R(6,1) = R^2*rho_R/muR * dUdr_R;
    if SurfBC==0
        u_R(2,1) = 0;
        u_R(4,1) = R/muR * tauTH;
    elseif SurfBC==1
        u_R(2,1) = vrTH;
        u_R(4,1) = 0;
    end
    % Check for consistency
%     u_R2 = PCR*u_C + PBR*du_B;
%     disp(double([u_R2(1) u_R(1)]))
%     disp(double([u_R2(2) u_R(2)]))
%     disp(double([u_R2(3) u_R(3)]))
%     disp(double([u_R2(4) u_R(4)]))
%     disp(double([u_R2(5) u_R(5)]))
%     disp(double([u_R2(6) u_R(6)]))
    
    
    for i=1:length(rp)
        r=rp(i)*R;
%         progressbar(i/length(rp))
        precise=1;
        [PBr,~,PRr]=ReturnPropagator6(Rvp(k),RCp,l,Visc,rp(i),precise);
        u_r = PRr * u_R;
        if rp(i)<Rvp(k)
            u_r = u_r - PBr*du_B;
        end
%         [PBR2,PCR2]=ReturnPropagator6(Rvp,RCp,l,Visc,rp(i));
%         u_r = PCR2*u_C;
%         if rp(i)>Rvp
%             u_r = u_r + PBR2*du_B;
%         end

        y1(i,ll,k) = u_r(1);
        y2(i,ll,k) = u_r(2);
        y3(i,ll,k) = u_r(3) * muR/r;
        y4(i,ll,k) = u_r(4) * muR/r;
        y5(i,ll,k) = u_r(5) * muR/rho_R/r;
        y6(i,ll,k) = u_r(6) * muR/rho_R/r^2;
    end
%         % Check for consistency        
%     [PBR2,PCR2]=ReturnPropagator6(Rvp(k),RCp,l,Visc,rp(end));
%     u_R = PCR2*u_C + PBR2*du_B;
%     disp(double([u_R(1) 0 y1(end)]))
%     disp(double([u_R(2) 0 y2(end)]))
%     disp(double([u_R(3) -R/muR*rho_m*gr*dH y3(end)]))
%     disp(double([u_R(4) tauTH y4(end)]))
%     disp(double([u_R(5) R*rho_m/muR * U_R y5(end)]))
%     disp(double([u_R(6) R^2*rho_m/muR * dUdr_R y6(end)]))

end
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


elseif strcmp(deg,'demo1')
    
    % Calculate state vector throughout the Moon
    R = 1737e3;
    CMB = 373e3;
    Rv = R-800e3;
    depths = 0:10e3:(R-Rv)+100e3;
    radii = fliplr(R*ones(size(depths))-depths);
    rp = radii/R
%     rp = vpa(1);
    deg = 10;
    Rvp = Rv/R;
    Visc='LunarVisc';
%     Visc='isoviscous';
    Te=0;
    SurfBC=1;
    [y1,y2,y3,y4,y5,y6]=dyn_flow6(deg,rp,Rvp,Visc,Te,'Moon',SurfBC);
    if length(y1)>1
        figure
        plot(radii,y1,'k')
        hold on
        plot(radii,y2,'r')
        title('Velocities')
        figure
        plot(radii,y3,'k')
        hold on
        plot(radii,y4,'r')
        title('Stresses')
%         figure
%         plot(radii,y5,'k')
%         figure
%         plot(radii,y6,'r')
    end

end

end

