function [lmcosi,dw]=xyz2plm(fthph,L,method,lat,lon,cnd)
% [lmcosi,dw]=XYZ2PLM(fthph,L,method,lat,lon,cnd)
%
% Forward real spherical harmonic transform in the 4pi normalized basis.
%
% Converts a spatially gridded field into spherical harmonics.
% For complete and regular spatial samplings [0 360 -90 90].
% If regularly spaced and complete, do not specify lat,lon.
% If not regularly spaced, fthph, lat and lon are column vectors.
%
% INPUT:
%
% fthph         Function defined on colatitude theta and longitude phi
% L             Maximum degree of the expansion (Nyquist checked)
% method        'gl'         By Gauss-Legendre integration (fast, inaccurate)
%               'simpson'    By Simpson integation (fast, inaccurate)
%               'irr'        By inversion (irregular samplings)
%               'im'         By inversion (fast, accurate, preferred)
%               'fib'        By Riemann sum on a Fibonacci grid (not done yet)
% lat           If not [90,-90], give latitudes explicitly, in degrees
% lon           If not [0,360], give longitudes explicitly, in degrees
% cnd           Eigenvalue tolerance in the irregular case
%
% OUTPUT:
%
% lmcosi        Matrix listing l,m,cosine and sine coefficients
% dw            Eigenvalue spectrum in the irregular case
%
% Method 'im' works with the orthogonality of the cosine and sine
% functions of the longitude instead of with the orthogonality of the
% Legendre polynomials. Longitudinal FFT is used in 'im', 'gl' and
% 'simpson', so longitudes must be equally spaced and periodic.
% Only method 'simpson' requires equidistant latitude spacing. For 'gl',
% the function is resampled to Gauss-Legende integration points. For
% 'im', can specify latitudes only, restriction that nlat>=(L+1).
%
% Note that the MEAN of the input data deviates from C(1), as sampled
% fields lose the orthogonality. The inversion approaches should recover
% the exact value of C(1), the true mean of the data, not the sample mean.
%
% See also PLM2XYZ, PLM2SPEC, PLOTPLM, etc.
%
% Last modified by fjsimons-at-alum.mit.edu, 11/05/2010

t0=clock;
warning off
defval('method','im')
defval('lon',[])
defval('lat',[])
defval('dw',[])
defval('cnd',[])

as=0;
% If no grid is specified, assumes equal spacing and complete grid
if isempty(lat) & isempty(lon)
  % Test if data is 2D, and periodic over longitude
  fthph=reduntest(fthph);
  chk=polestest(fthph);
%   if chk==0
%       chk=0;  
%   end
  % Make a complete grid
  nlon=size(fthph,2);
  nlat=size(fthph,1);
  % Nyquist wavelength
  Lnyq=min([ceil((nlon-1)/2) nlat-1]);
  % Colatitude and its increment
  theta=linspace(0,pi,nlat);
  as=1; % Equally spaced
  % Calculate latitude/longitude sampling interval; no wrap-around left
  dtheta=pi/(nlat-1);
  dphi=2*pi/nlon;
  switch method 
    % Even without lat/lon can still choose the full inversion method
    % without Fourier transformation
   case 'irr'
     [LON,LAT]=meshgrid(linspace(0,2*pi*(1-1/nlon),nlon),...
			linspace(pi/2,-pi/2,nlat));
     lat=LAT(:); lon=LON(:); fthph=fthph(:);
     theta=pi/2-lat;
     clear LON LAT
  end
elseif isempty(lon)
  % If only latitudes are specified; make equal spacing longitude grid
  % Latitudes can be unequally spaced for 'im', 'irr' and 'gl'.
  fthph=reduntest(fthph);
  theta=(90-lat)*pi/180;
  dtheta=(lat(1)-lat(2))*pi/180;
  nlat=length(lat);
  nlon=size(fthph,2);
  dphi=2*pi/nlon;
  Lnyq=min([ceil((nlon-1)/2) ceil(pi/dtheta)]);
else
  % Irregularly sampled data
  fthph=fthph(:);
  theta=(90-lat)*pi/180;
  lat=lat(:)*pi/180;
  lon=lon(:)*pi/180;
  nlon=length(lon);
  nlat=length(lat);
  % Nyquist wavelength
  adi=[abs(diff(sort(lon))) ; abs(diff(sort(lat)))];
  Lnyq=ceil(pi/min(adi(~~adi)));
  method='irr';
end

% Decide on the Nyquist frequency
defval('L',Lnyq);
% Never use Libbrecht algorithm... found out it wasn't that good
defval('libb',0)
%disp(sprintf('Lnyq= %i ; expansion out to degree L= %i',Lnyq,L))

if L>Lnyq | nlat<(L+1)
  error('XYZ2PLM: Function undersampled. Aliasing will occur.')
end

% Make cosine and sine matrices
[m,l,mz]=addmon(L);
lmcosi=[l m zeros(length(l),2)];

% Define evaluation points
switch method
 case 'gl'
  % Highest degree of integrand will always be 2*L
  [w,x]=gausslegendrecof(2*L,[],[-1 1]);
  % Function interpolated at Gauss-Legendre latitudes; 2D no help
  fthph=interp1(theta,fthph,acos(x),'spline');
 case {'irr','simpson','im'}
  % Where else to evaluate the Legendre polynomials
  x=cos(theta);
 otherwise
  error('Specify valid method')
end

fnpl=sprintf('%s/LSSM-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'LEGENDRE'),L,length(x));
 
if exist(fnpl,'file')==2 & as==1
  load(fnpl)
else  
  % Evaluate Legendre polynomials at selected points
  Plm=repmat(NaN,length(x),addmup(L));
  if L>200
    h=waitbar(0,'Evaluating all Legendre polynomials');
  end
  in1=0;
  in2=1;
  for l=0:L
    if libb==0
      Plm(:,in1+1:in2)=(legendre(l,x(:)','sch')*sqrt(2*l+1))';
    else
      Plm(:,in1+1:in2)=(libbrecht(l,x(:)','sch')*sqrt(2*l+1))';
    end
    in1=in2;
    in2=in1+l+2;
    if L>200
      waitbar((l+1)/(L+1),h)
    end
  end
  if L>200
    delete(h)
  end
  if as==1
    save(fnpl,'Plm','-v7.3')
  end
end

switch method
 case {'irr'}
  Plm=[Plm.*cos(lon(:)*m(:)') Plm.*sin(lon(:)*m(:)')];
  % Add these into the sensitivity matrix
  [C,merr,mcov,chi2,L2err,rnk,dw]=datafit(Plm,fthph,[],[],cnd);
  lmcosi(:,3)=C(1:end/2);
  lmcosi(:,4)=C(end/2+1:end);
 case {'im','gl','simpson'}
  % Perhaps demean the data for Fourier transform
  defval('dem',0)
  if dem
    meanm=mean(fthph,2);
    fthph=fthph-repmat(meanm,1,nlon);
  end
  
  % Calculate integration over phi by the fast Fourier
  % transform. Integration of real input field with respect to the second
  % dimension of r, at  wavenumber m, thus at constant latitude. You get
  % as many wavenumbers m as there are longitudes; only use to L. With
  % Matlab's FFT, need to multiply by sampling interval. 
  gfft=dphi*fft(fthph,nlon,2);

  if dem
    % Add the azimuthal mean back in there
    gfft(:,1)=2*pi*meanm;
  end

  % Note these things are only half unique - the maximum m is nlon/2
  % But no Nyquist theory exists for the Legendre transform...
  a=real(gfft);
  b=-imag(gfft);
  in1=0;
  in2=1;
 otherwise
  error('Specify valid method')
end

switch method
 case 'im'
  % Loop over the orders. This speeds it up versus 'irr'
  for ord=0:L
    a(:,1)=a(:,1)/2;
    b(:,1)=b(:,1)/2;
    Pm=Plm(:,mz(ord+1:end)+ord)*pi;
    [lmcosi(mz(ord+1:end)+ord,3)]=datafit(Pm,a(:,ord+1),[],[],cnd);
    [lmcosi(mz(ord+1:end)+ord,4)]=datafit(Pm,b(:,ord+1),[],[],cnd);
  end
 case 'simpson'
  % Loop over the degrees. Could go up to l=nlon if you want
  for l=0:L,
    % Integrate over theta using Simpson's rule
    clm=simpson(theta,...
		repmat(sin(theta(:)),1,l+1).*a(:,1:l+1).*Plm(:,in1+1:in2));
    slm=simpson(theta,...
		repmat(sin(theta(:)),1,l+1).*b(:,1:l+1).*Plm(:,in1+1: ...
						  in2));
    in1=in2;
    in2=in1+l+2;
    % And stick it in a matrix [l m Ccos Csin]
    lmcosi(addmup(l-1)+1:addmup(l),3)=clm(:)/4/pi;
    lmcosi(addmup(l-1)+1:addmup(l),4)=slm(:)/4/pi;
  end
 case 'gl'
  % Loop over the degrees. Could go up to l=nlon if you want
  for l=0:L,
    % Integrate over theta using Gauss-Legendre integration
    clm=sum(a(:,1:l+1).*(diag(w)*Plm(:,in1+1:in2)));
    slm=sum(b(:,1:l+1).*(diag(w)*Plm(:,in1+1:in2)));
    in1=in2;
    in2=in1+l+2;
    % And stick it in a matrix [l m Ccos Csin]
    lmcosi(addmup(l-1)+1:addmup(l),3)=clm(:)/4/pi;
    lmcosi(addmup(l-1)+1:addmup(l),4)=slm(:)/4/pi;
  end
  rnk=[]; 
end

% Get rid of machine precision error
lmcosi(abs(lmcosi(:,3))<eps,3)=0;
lmcosi(abs(lmcosi(:,4))<eps,4)=0;

%disp(sprintf('XYZ2PLM (Analysis)  took %8.4f s',etime(clock,t0)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grd=reduntest(grd)
% Tests if last longitude repeats last (0,360)
% and removes last data column
if sum(abs( (grd(:,1)-grd(:,end))/mean(grd(:,1)) )) >= size(grd,2)*eps*1e6
  disp(sprintf('Data violate wrap-around by %8.4e',...
		  sum(abs(grd(:,1)-grd(:,end)))))
end
grd=grd(:,1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chk=polestest(grd)
% Tests if poles (-90,90) are identical over longitudes 
chk=1;
var1=var(grd(1,:))/mean(grd(1,:));
var2=var(grd(end,:))/mean(grd(1,:));
if var1>eps*1e6 || var2>eps*1e6
    chk=0;
    disp(sprintf('Poles violated by %8.4e and %8.4e',var1,var2))
end


