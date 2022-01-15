function [lmcosi,finite,H0]=...
    Topo2Grav(hlmcosi,drho,r,rp,GM,RH,wat,npower,lmax,lmaxh)
% [lmcosi,finite]=...
%    Topo2Grav(hlmcosi,drho,r,rp,GM,R,wat,npower,lmax,lmaxh)
%
% This function calculates gravity at radius r from finite-amplitude 
% interface relief hlmcosi at radius rp.
%   
% OUTPUT:
% lmcosi     Complete gravity coefficients
% finite     Only the finite amplitude portion
% H0         Adjustment to rp (useful if input is a grid)
%
% INPUT:
% hlmcosi    inputed interface relief [m], either lmcosi array or lat-lon
%            gridded data
% drho       density contrast, either a scalar, lmcosi, or a map with the
%            same resolution as the topography map (as dictated by lmaxh)
% r          Measurement radius [m]
% rp         Interface radius [m]   NOTE: if hlmcosi(1,3) contains the
%            interface radius, this input is redundant.
% GM         Gravitational constant times mass [m^3/s^2] (alternatively, 
%            planet name)
% RH         Planetary radius
% wat        0 geoid anomaly [m]
%            1 gravitational potential anomaly [J/kg] - Default
%            2 free-air gravity anomaly [m/s^2]
%            3 gravity disturbance [m/s^2]
%            4 Clm coefficients [1]
% npower     max power of topography for finite amplitude calculation
% lmax       maximum degree of the output
% lmaxh      maximum degree of the topographic power calculation (used for
%            calculating dres)
% 
% Last modified 10/10/16 by pjames at alum dot mit dot edu

Grav=6.67408*10^-11;
defval('drho',500)
defval('wat',1) % Gravitational potential
defval('r',1738e3); % Moon
% defval('rp',r);
defval('GM','Moon');
defval('npower',6);
defval('up',1);
defval('lmax',hlmcosi(end,1)) % Assumes hlmcosi is an lmcosi matrix!
defval('lmaxh',2*lmax)

if ischar(GM)
    if strcmp(GM,'Moon')
        GM=4.9028e+12;
        RH=1738e3;
    elseif strcmp(GM,'Mercury')
        GM=2.2032e13;
        RH=2440e3;
    elseif strcmp(GM,'Mars')
        GM=4.282837e13;
        RH=3390e3;
    elseif strcmp(GM,'Venus')
        GM=3.24859e13;
        RH=6052e3;
    elseif strcmp(GM,'Earth')
        GM=3.986004418e13;
        RH=6371e3;
    end
end
M=GM/Grav;
% gacc=GM/R^2;
H0=0; % The degree zero term of topography
rhoH0=0; 
ChangeRadius=1;
if r>=rp
    up=1; 
else
    up=0; 
end

if ~ischar(hlmcosi)
    
[~,~,~,Clm]=addmon(lmax);
[~,~,~,Clm_finite]=addmon(lmax);

if size(hlmcosi,2)>6 % Assume the input is a lat-lon grid
    hmap=hlmcosi;
    dres=180/(size(hmap,1)-1);
    hlmcosi=xyz2plm(hmap,lmax); % If the inputted lmax is too high, this 
                                % will produce an error
    % The lmcosi array calculated from hmap might have a 0,0 term. 
    if hlmcosi(1,3)~=0
        if ChangeRadius==1
            if rp~=hlmcosi(1,3) 
                if hlmcosi(1,3)/rp>0.5
                    % Probably a new radius
                    rp=hlmcosi(1,3);
                    hmap = hmap - rp*ones(size(hmap));
                else
                    % Probably an adjustment to the radius
                    H0 = hlmcosi(1,3);
                    rp = rp + H0;
                    hmap = hmap - H0*ones(size(hmap));
                    disp('Adjustment:')
                    disp(hlmcosi(1,3))
                end
                disp('New Rp:')
                disp(rp)
            else
                % Leave the radius the same
                disp('Radius unchanged:')
                disp(rp)
                hmap = hmap - rp*ones(size(hmap));
            end
        end
        hlmcosi(1,3) = 0;
    end
else % Assume the input is an lmcosi array
    if hlmcosi(1,3)~=0
        if ChangeRadius==1
            if rp~=hlmcosi(1,3)
                if hlmcosi(1,3)/rp>0.5
                    % Probably a new radius
                    rp=hlmcosi(1,3);
                else
                    % Probably an adjustment to the radius
                    rp=rp+hlmcosi(1,3);
                end
                disp('New value for Rp:')
                disp(rp)
            end
        end
        H0=hlmcosi(1,3);
        hlmcosi(1,3)=0;
    end
    hlmcosi=hlmcosi(1:addmup(lmax),1:4);
    if size(drho,2)>6 % Check if drho is a map
        % Find dres to match drho
        dres=180/(size(drho,1)-1);
    else
        dres=(180/sqrt(lmax*(lmax+1))) * lmax/lmaxh;
    end
    hmap=plm2xyz(hlmcosi,dres);
end

% If rho varies spatially, calculate maps for topography and density
if size(drho,2)>1
    if size(drho,2)>6
        % Assume the input is a map
        if size(drho)~=size(hmap)
            warning('density map doesn''t match topography map')
        end
        rhomap=drho;
        rhohlmcosi=xyz2plm(rhomap.*hmap,lmax);
    else
        % Assume the input is an lmcosi array
%         rholmcosi=drho;
        rhomap=plm2xyz(drho,dres);
        rhohlmcosi=xyz2plm(rhomap.*hmap,lmax);
    end
    rhoh0=rhohlmcosi(1,3);
else
    % drho is a constant density contrast
    rhomap=drho*ones(size(hmap));
    rhoh0=0;
end

% First order gravity (without finite amplitude effects):
for i=2:size(hlmcosi,1)
    l=hlmcosi(i,1);
    if size(drho,2)>1
        Clm(i,3:4) = 4*pi*rp^2/(M*(2*l+1)) * rhohlmcosi(i,3:4);
    else
        Clm(i,3:4) = 4*pi*rp^2/(M*(2*l+1)) * drho*hlmcosi(i,3:4);
    end
end
% Calculate finite amplitude portion of gravity
if npower>1
%     hmap=plm2xyz(hlmcosi,dres);
    for n=2:npower
%         if size(drho,2)>1
            rhohpower=xyz2plm(hmap.^n.*rhomap,lmax);
%         else
%             rhohpower=xyz2plm(hmap.^n*drho,lmax);
%         end
        for i=2:size(hlmcosi,1)
            l=hlmcosi(i,1);
            Prod=1;
            for j=1:n
                if up
                    Prod=Prod*(l+4-j)/j/rp;
                else
                    Prod=Prod*(l+j-3)/j/rp;
                end
            end
            if up
                Prod=Prod/(l+3);
            else
                if l==2
                    % Singularity, re-do the loop over powers
                    % j==1:
                    Prod=1/rp;
                    if n>1
                        for j=2:n
                            Prod=Prod*(l+j-3)/j/rp;
                        end
                    end
                else
                Prod=Prod/(l-2);
                end
            end

            Clm_finite(i,3:4) = Clm_finite(i,3:4) + ...
                4*pi*rp^3/(M*(2*l+1)) * rhohpower(i,3:4) * Prod;
        end
    end
end

% Add finite amplitude portion to the gravity coefficients
Clm(:,3:4) = Clm(:,3:4) + Clm_finite(:,3:4);

for i=2:size(Clm,1)
    l=Clm(i,1);
    if up==1
        % Upward continuation
        Clm(i,3:4)=Clm(i,3:4)*(rp/r)^l;
        Clm_finite(i,3:4)=Clm_finite(i,3:4)*(rp/r)^l;
    elseif up==0
        % Downward continuation
        Clm(i,3:4)=Clm(i,3:4)*(r/rp)^(l+1);
        Clm_finite(i,3:4)=Clm_finite(i,3:4)*(r/rp)^(l+1);
    end
end

% Convert to desired form of gravity (comparable to plm2pot)
lmcosi=Clm; finite=Clm_finite;
for i=2:size(lmcosi,1)
    l=lmcosi(i,1);
    if wat==0 % geoid
        lmcosi(i,3:4)=Clm(i,3:4)*r;
        finite(i,3:4)=Clm_finite(i,3:4)*r;
    elseif wat==1 % Gravitational potential
        lmcosi(i,3:4)=Clm(i,3:4)*GM/r;
        finite(i,3:4)=Clm_finite(i,3:4)*GM/r;
    elseif wat==2 % Free-air gravity
        lmcosi(i,3:4)=Clm(i,3:4)*(l-1)*GM/r^2;
        finite(i,3:4)=Clm_finite(i,3:4)*(l-1)*GM/r^2;
    elseif wat==3 % Gravity disturbance
        lmcosi(i,3:4)=Clm(i,3:4)*(l-1)*GM/r^2;
        finite(i,3:4)=Clm_finite(i,3:4)*(l-1)*GM/r^2;
    elseif wat==4 % coefficients
        lmcosi(i,3:4)=Clm(i,3:4);
        finite(i,3:4)=Clm_finite(i,3:4);
    else
        warning('ON','Input not understood')
    end
end

Clm_nofinite = Clm; Clm_nofinite(:,3:4) = Clm(:,3:4)-Clm_finite(:,3:4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(hlmcosi,'demo1')
    
    lmax=360;

    R=1738e3;

    H = load('SHAPE/LOLA1439v2.sh');
    H = H(1:addmup(lmax),1:4);
    a = H(1,3);
    H(1,3)=a-R;
    % Greg's gravity from topography, referenced to 1738 km
    GFTlm = load('GRAV/FROMTOPO/LOLA1500_PA_MWRP_1000KGM3_N14_SHA.TAB');
    GFTlm = GFTlm(1:addmup(lmax),1:4); GFTlm(:,3:4) = 2.56*GFTlm(:,3:4);
    [Clm,GM,R] = read_TAB('GRAV/JGGRAIL_900C11A_SHA.TAB');
    Clm = Clm(1:addmup(lmax),1:4);
    [Blm,~,~] = read_TAB('GRAV/JGGRAIL_900C11A_BOUGUER_SHA.TAB');
    Blm = vertcat([0 0 0 0; 1 0 0 0; 1 1 0 0; 2 0 0 0],Blm);
    Blm = Blm(1:addmup(lmax),1:4);
    N = plm2pot(Clm,R,GM,R,0,'spherical');
    %
%     plotmap(H); title('Topography'); colorbar
%     plotmap(Clm); title('Coefficients'); colorbar
%     plotmap(GFTlm); title('Gravity from Topography'); colorbar
%     plotmap(N); title('Geoid'); colorbar

    % Reproduce gravity from topography

    drho=2560;
%     Rp=R;
    wat=4; % Clm coefficients
    npower=6;
    a
    R
%     if R>=a; up=1; else up=0; end
    [GFTlm_check,GFTlm_finite] = Topo2Grav(H,drho,R,a,GM,R,wat,npower);

    GDFTcheck=plm2pot(GFTlm_check,R,GM,R,3,'spherical');
    GDFTfinite=plm2pot(GFTlm_finite,R,GM,R,3,'spherical');
    GDFT=plm2pot(GFTlm,R,GM,R,3,'spherical');

%     plotmap(GFTlm_check); title('Calculated Gravity from Topo'); colorbar
%     plotmap(GFTlm_finite); title('Finite amplitude coefficients'); colorbar
%     plotmap(GDFTfinite); title('Finite amplitude gravity'); colorbar
    plotmap(GDFTcheck); title('Calculated gravity from topo'); colorbar
    plotmap(GDFT); title('Actual gravity from topo'); colorbar
    plotmap([GDFT(:,1:2) GDFT(:,3:4)-GDFTcheck(:,3:4)])
    
elseif strcmp(hlmcosi,'demo2')    
    lmax=180;

    R=1738e3;

    H = load('SHAPE/LOLA1439v2.sh');
    H = H(1:addmup(lmax),1:4);
    a = H(1,3);
    H(1,3)=a-R;
    GFTlm = load('GRAV/FROMTOPO/LOLA1500_PA_MWRP_1000KGM3_N14_SHA.TAB');
    GFTlm = GFTlm(1:addmup(lmax),1:4); GFTlm(:,3:4) = 2.56*GFTlm(:,3:4);
    [Clm,GM,R] = read_TAB('GRAV/JGGRAIL_900C11A_SHA.TAB');
    Clm = Clm(1:addmup(lmax),1:4);
    [Blm,~,~] = read_TAB('GRAV/JGGRAIL_900C11A_BOUGUER_SHA.TAB');
    Blm = vertcat([0 0 0 0; 1 0 0 0; 1 1 0 0; 2 0 0 0],Blm);
    Blm = Blm(1:addmup(lmax),1:4);

    drho=2560;
    Rp=R;
    wat=4; % Clm coefficients
    npower=6;
    up=1;

    [GFTlm_check,GFTlm_finite] = Topo2Grav(H,drho,R,Rp,GM,R,wat,npower);

    Clm_residual=Clm;
    Clm_residual(:,3:4)=Clm_residual(:,3:4)-GFTlm_check(:,3:4);

    Gresidual=plm2pot(Clm_residual,R,GM,R,3,'spherical');
    Gboug=plm2pot(Blm,R,GM,R,3,'spherical');
    Gresidual2=Gresidual; Gresidual2(1:4,3:4)=zeros(4,2);

    plotmap(Blm); title('Bouguer gravity coefficients'); colorbar
    plotmap(Gresidual); title('Calculated Bouguer gravity'); colorbar
    plotmap(Gresidual2); title('Calculated Bouguer w/o low degrees'); colorbar
    plotmap(Gboug); title('Bouguer gravity'); colorbar

    % It works, but Gresidual includes degree 1 terms and (2,0), whereas Gboug
    % does not...

elseif strcmp(hlmcosi,'demo3')
    % Compare to Anton's code, and compare to Greg's gravity from LOLA
    
    lmax=90;
    Grav=6.67408*10^-11;
    R=1738e3;
    npower=6;
%     lmaxh=npower*lmax;
    lmaxh=min([npower 3])*lmax;

    H = load('SHAPE/LOLA1439v2.sh');
    H = H(1:addmup(lmax),1:4);
    a = H(1,3);
    drho=2560;
    Hmap=plm2xyz(H,180/lmaxh);
%     plotmap(Hmap); title('Shape')

    % Greg's gravity from topography, referenced to 1738 km
    GFTlm = load('GRAV/FROMTOPO/LOLA1500_PA_MWRP_1000KGM3_N14_SHA.TAB');
    GFTlm = GFTlm(1:addmup(lmax),1:4); 
    GFTlm(:,3:4) = drho/1000*GFTlm(:,3:4);

    wat=4; % Coefficients        
%     [lmcosi_topo,~] = xyz2plm(Hmap,nmaxgt,'im',[],[],[]);

    GFTlm_Peter = Topo2Grav(H,drho,R,[],'Moon',[],wat,npower,lmax,lmaxh);
%     GFTlm_Peter = Topo2Grav(lmcosi_topo,drho,R,[],'Moon',wat,npower,lmax,lmax*3);



    
    lmcosi_grav = Topo2Grav_Anton(Hmap,R,lmaxh,lmax,npower);
    
    rhobar=GM/Grav / (4/3*pi*a^3);
    lmcosi_grav(:,3:4) = lmcosi_grav(:,3:4) * drho/rhobar;
    
    lmcosi_grav(1,3)=0;
    
    wat=1; % Gravity disturbance
    GD_Anton=plm2pot(lmcosi_grav,R,GM,R,wat,'spherical');
    GD_Peter=plm2pot(GFTlm_Peter,R,GM,R,wat,'spherical');
    GD_Greg=plm2pot(GFTlm,R,GM,R,wat,'spherical');
    Sff_Anton=calculate_power(GD_Anton);
    Sff_Peter=calculate_power(GD_Peter);
    Sff_Greg=calculate_power(GD_Greg);
    
    figure
    semilogy(Sff_Anton,'r')
    hold on
    semilogy(Sff_Peter,'k')
    semilogy(Sff_Greg,'g')
    figure
    plot(Sff_Anton./Sff_Peter)
    % ylim([1.705 1.71])
    figure
    plot(Sff_Peter./Sff_Greg)

    Greg_map=plm2xyz(GD_Greg);
    Peter_map=plm2xyz(GD_Peter);
    Anton_map=plm2xyz(GD_Anton);
    sf=mean(mean(Peter_map./Anton_map));
    
    plotmap(Greg_map); title('Greg''s gravity from topo'); colorbar
    plotmap(Peter_map); title('Peter''s gravity from topo'); colorbar
    plotmap(Anton_map); title('Anton''s gravity from topo'); colorbar
    plotmap(Greg_map-Peter_map); title('Greg-Peter'); colorbar
    plotmap(Greg_map-Anton_map); title('Greg-Anton'); colorbar
    plotmap(Peter_map-Anton_map); title('Peter-Anton'); colorbar
    
else warning('ON','Input not understood')

end

end


