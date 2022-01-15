function [DidItConverge,CF]=CheckConvergence(model,oldmodel,criterion)

% Check to see if a spherical harmonic model is converging
% INPUT:
% model         -   Current model, lmcosi coefficients or a map
% oldmodel      -   Previous model
% criterion     -   Threshold factor by which the iterations are changing
% OUTPUT:
% DidItConverge -   1: It did!
%                   0: Not yet...
%                   2: It may be blowing up

defval('criterion',1e-6)
warning ON

if size(model,2)<7 && size(model,2)>3 % Spherical harmonic coefficients
    model_delta=model; 
    model_delta(:,3:4)=model(:,3:4)-oldmodel(:,3:4);
    
    lmax=max(model(:,1));
    Smm=zeros(lmax,1);
    Sdd=zeros(lmax,1);
    for i=1:length(model(:,1))
        l=model(i,1);
        if l>0
            Smm(l)=Smm(l)+model(i,3)^2+model(i,4)^2;         
            Sdd(l)=Sdd(l)+model_delta(i,3)^2+model_delta(i,4)^2;            
        end
    end


    CF=max(Sdd./Smm);
    if CF<criterion
        DidItConverge=1;
        disp(['Convergence factor = ' num2str(CF)])
    else
        DidItConverge=0;
    end
    
    % What could go wrong?
    % ...High-degree power could massively exceed low-degree power
    lmax=length(Smm);
    if max(Smm(round(lmax/2):lmax))>mean(Smm(1:5)*1000)
        warning('Elevated power at high degrees')
        DidItConverge=2;
    end        

else
    warning('Please input spherical harmonic coefficients')
end



end