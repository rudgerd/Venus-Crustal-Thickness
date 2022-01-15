function lmcosi=Density2Grav(rlmcosi,r,R1,R2,wat,gacc)

% This function calculates potential perturbation from lateral variations
% in density (see Wieczorek at al., 2013, SOM, section 5)
%   
% OUTPUT:
% lmcosi     complete gravity coefficients
%
% INPUT:
% rlmcosi    inputed crustal density [m]
% r          measurement radius [m]
% R1         upper radius of density layer [m]
% R2         lower radius of density layer [m]
% wat        0 - Geoid coefficients [m]
%            1 - Potential coefficients [m^2/s^2] - [Default]
%            2 - Free air gravity coefficients [m/s^2]
%            3 - Gravity Disturbance [m/s^2]
%            4 - Coefficients [1]
% gacc       gravitational acceleration at radius r [m/s^2]
% 
% Last modified 4/29/16 by pjames at alum dot mit dot edu

Grav=6.67300*10^-11;
defval('drho',500)
defval('wat',1)
defval('r',2440e3); % Mercury
defval('R1',r);
defval('R2',R1-40e3);
defval('gacc',3.7);
defval('npower',6);
defval('up',1);


if ~ischar(rlmcosi)

if R1<R2; warning('ON','Radius R1 must be larger than radius R2'); end
if r<R1 && r>R2; warning('ON','Sampling radius cannot be between R2 and R1'); end
if length(rlmcosi)<2; warning('ON','Input rlmcosi is a scalar'); end

lmax=rlmcosi(end,1);
[~,~,~,Ulmcosi]=addmon(lmax);
% Sum over powers of topography
for i=2:size(rlmcosi,1);
    l=rlmcosi(i,1);
    if r>=R1 % Upward continuation
        Ulmcosi(i,3:4) = 4*pi*Grav*r^2/(2*l+1)/(l+3) * ...
            ((R1/r)^(l+3)-(R2/r)^(l+3)) * rlmcosi(i,3:4);
    else % Downward continuation
        if l==2
            Ulmcosi(i,3:4) = 4*pi*Grav*r^2/(2*l+1) * ...
                log(R1/R2) * rlmcosi(i,3:4);
        else
            Ulmcosi(i,3:4) = 4*pi*Grav*r^2/(2*l+1)/(l-2) * ...
                ((r/R2)^(l-2)-(r/R1)^(l-2)) * rlmcosi(i,3:4);
        end
    end
    
end

% Convert to desired form of gravity (comparable to plm2pot)
lmcosi=Ulmcosi;
for i=2:size(lmcosi,1);
    if wat==0 % geoid
        lmcosi(i,3:4)=Ulmcosi(i,3:4)/gacc;
    elseif wat==1 % Gravitational potential
    elseif wat==2 % Free-air gravity
        lmcosi(i,3:4)=Ulmcosi(i,3:4)*(l-1)/r;
    elseif wat==3 % Gravity disturbance
        lmcosi(i,3:4)=Ulmcosi(i,3:4)*(l-1)/r;
    elseif wat==4 % coefficients
        lmcosi(i,3:4)=Ulmcosi(i,3:4)/r/gacc;
    else
        warning('ON','Input not understood')
    end
end

elseif strcmp(rlmcosi,'demo1')
    
    GrainDensity=load('DATA/grain_density_310.sh');
    porosity=.12;
    rlmcosi=GrainDensity; rlmcosi(:,3:4)=rlmcosi(:,3:4)*(1-porosity);
    lmax=GrainDensity(end,1);

    R=1738e3;
    Rp=R-43e3;
    gacc=1.6;
    wat=0; % Geoid
    
    Nlmcosi=Density2Grav(rlmcosi,R,R,Rp,wat,gacc);
%     H = load('SHAPE/LOLA1439v2.sh');
%     H = H(1:addmup(lmax),1:4);
%     a = H(1,3);
%     H(1,3)=a-R;
%     GFTlm = load('GRAV/FROMTOPO/LOLA1500_PA_MWRP_1000KGM3_N14_SHA.TAB');
%     GFTlm = GFTlm(1:addmup(lmax),1:4); GFTlm(:,3:4) = 2.56*GFTlm(:,3:4);
%     [Clm,GM,R] = read_TAB('GRAV/JGGRAIL_900C11A_SHA.TAB');
%     Clm = Clm(1:addmup(lmax),1:4);
%     [Blm,~,~] = read_TAB('GRAV/JGGRAIL_900C11A_BOUGUER_SHA.TAB');
%     Blm = vertcat([0 0 0 0; 1 0 0 0; 1 1 0 0; 2 0 0 0],Blm);
%     Blm = Blm(1:addmup(lmax),1:4);
%     N = plm2pot(Clm,R,GM,R,0,'spherical');
    %
    plotmap(GrainDensity); title('Density'); colorbar
    plotmap(rlmcosi); title('Density with porosity'); colorbar
    plotmap(Nlmcosi); title('Geoid'); colorbar

    Nlmcosi2=Density2Grav(rlmcosi,Rp,R,Rp,wat,gacc);
    plotmap(Nlmcosi2); title('Geoid'); colorbar
%%
    plotmap(Nlmcosi2(addmup(1)-1:addmup(2),:),[],1); 
    title('Geoid degrees 1 and 2');
    plotmap(rlmcosi(addmup(1)-1:addmup(2),:),[],1); 
    title('Density degrees 1 and 2');
    Nlmcosi3=Nlmcosi2; Nlmcosi3(4:6,3:4)=[0,0;0,0;0,0];
    plotmap(Nlmcosi3); title('Geoid without degrees 1 or 2');

else warning('ON','Input not understood')

end

end
