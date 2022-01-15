function [PBr,PCr,PRr]=ReturnPropagator6(B,C,l,MantleProperties,rp,precise)
%  Mantle Properties: a nx3 matrix with non-dimensional radii in the first
%  column, viscosities in the second column, and densities (if desired) in
%  the third column.  Radii in the first column should increase
%
%  rp - The radius to which the state vector is propagated (if propagating
%  to the planetary surface, this equals zero.

defval('rp',1)
defval('precise',0)
R=1; % Planetary radius, normalized by the planetary radius

% Check the input arguments
if max(MantleProperties(:,1)) > 1
   warning('Interface radii should be normalized');
elseif norm(MantleProperties(:,1)-sort(MantleProperties(:,1))) > 1E-8
   warning('Interface radii in MantleProperties are not correctly sorted');
   disp(MantleProperties)
elseif MantleProperties(1,1) < C
   warning('Mantle Properties are specified deeper than the lower BC');
   disp(MantleProperties(1,1))
   disp(C)
elseif MantleProperties(size(MantleProperties,1)) ~= 1
   warning('Properties should be specified all the way to the surface');
elseif B <= C
   warning('Mantle load should be above the lower BC')
end

if precise==0
    % Outputs that propagate upward
    PCr=eye(6); if B<rp; PBr=eye(6); end
    for i=1:size(MantleProperties,1)
        rlayer=MantleProperties(i,1);
        if i==1
            rlayer2=C;
        else
            rlayer2=MantleProperties(i-1,1);
        end

        if rlayer2>=rp % Stop propagating
            % Do nothing
        else
            if rlayer>rp
                % Propagate once more
                rlayer=rp;
            end
            mustar=MantleProperties(i,2);
            rhostar=MantleProperties(i,3);
            % Variable order:
            % vr
            % vtheta
            % tau_rr
            % tau_rtheta
            % U
            % dUdr

            L=l*(l+1);
            A = [-2, L, 0, 0, 0, 0;...
                -1, 1, 0, 1/mustar, 0, 0;...
                12*mustar, -6*L*mustar, 1, L, 0, -rhostar;...
                -6*mustar, 2*(2*L-1)*mustar, -1, -2, -rhostar, 0;...
                0, 0, 0, 0, 1, 1;...
                0, 0, 0, 0, L, 0];
            % Make PCr
            PCr=expm((log(rlayer)-log(rlayer2))*A)*PCr;
            if B<rp
                % Make PBr
                if B<rlayer
                    if B>rlayer2
                        PBr=expm((log(rlayer)-log(B))*A)*PBr;
                    else
                        PBr=expm((log(rlayer)-log(rlayer2))*A)*PBr;
                    end
                end
            end
        end
    end
    
    % Outputs that propagate downward:
    PRr=eye(6); if B>=rp; PBr=eye(6); end
    for i=1:size(MantleProperties,1)
        MP_down=flipud(MantleProperties);
        rlayer2=MP_down(i,1);
        if i==size(MP_down,1)
            rlayer=C;
        else
            rlayer=MP_down(i+1,1);
        end

        if rlayer2<=rp % Stop propagating
            % Do nothing
        else
            if rlayer<rp
                % Propagate once more
                rlayer=rp;
            end
            mustar=MP_down(i,2);
            rhostar=MP_down(i,3);
            L=l*(l+1);
            A = [-2, L, 0, 0, 0, 0;...
                -1, 1, 0, 1/mustar, 0, 0;...
                12*mustar, -6*L*mustar, 1, L, 0, -rhostar;...
                -6*mustar, 2*(2*L-1)*mustar, -1, -2, -rhostar, 0;...
                0, 0, 0, 0, 1, 1;...
                0, 0, 0, 0, L, 0];

            % Make PRr
            PRr=expm((log(rlayer)-log(rlayer2))*A)*PRr;
            if B>=rp
                % Make PBr
                if B<rlayer2
                    if B>rlayer
                        PBr=expm((log(rlayer)-log(B))*A)*PBr;
                    else
                        PBr=expm((log(rlayer)-log(rlayer2))*A)*PBr;
                    end
                end
            end
        end
    end
    
elseif precise==1
    % Outputs that propagate upward
    PCr=eye(6); if B<rp; PBr=eye(6); end
    for i=1:size(MantleProperties,1)
        rlayer=vpa(MantleProperties(i,1));
        if i==1
            rlayer2=vpa(C);
        else
            rlayer2=vpa(MantleProperties(i-1,1));
        end

        if rlayer2>=rp % Stop propagating
            % Do nothing
        else
            if rlayer>rp
                % Propagate once more
                rlayer=rp;
            end
            mustar=vpa(MantleProperties(i,2));
            rhostar=vpa(MantleProperties(i,3));
            % Variable order:
            % vr
            % vtheta
            % tau_rr
            % tau_rtheta
            % U
            % dUdr

            L=l*(l+1);
            A = vpa([-2, L, 0, 0, 0, 0;...
                -1, 1, 0, 1/mustar, 0, 0;...
                12*mustar, -6*L*mustar, 1, L, 0, -rhostar;...
                -6*mustar, 2*(2*L-1)*mustar, -1, -2, -rhostar, 0;...
                0, 0, 0, 0, 1, 1;...
                0, 0, 0, 0, L, 0]);
            % Make PCr
            PCr=expm((log(rlayer)-log(rlayer2))*A)*PCr;
            if B<rp
                % Make PBr
                if B<rlayer
                    if B>rlayer2
                        PBr=expm((log(rlayer)-log(B))*A)*PBr;
                    else
                        PBr=expm((log(rlayer)-log(rlayer2))*A)*PBr;
                    end
                end
            end
        end
    end
    
    % Outputs that propagate downward:
    PRr=eye(6); if B>=rp; PBr=eye(6); end
    for i=1:size(MantleProperties,1)
        MP_down=flipud(MantleProperties);
        rlayer2=vpa(MP_down(i,1));
        if i==size(MP_down,1)
            rlayer=vpa(C);
        else
            rlayer=vpa(MP_down(i+1,1));
        end

        if rlayer2<=rp % Stop propagating
            % Do nothing
        else
            if rlayer<rp
                % Propagate once more
                rlayer=rp;
            end
            mustar=vpa(MP_down(i,2));
            rhostar=vpa(MP_down(i,3));
            L=l*(l+1);
            A = vpa([-2, L, 0, 0, 0, 0;...
                -1, 1, 0, 1/mustar, 0, 0;...
                12*mustar, -6*L*mustar, 1, L, 0, -rhostar;...
                -6*mustar, 2*(2*L-1)*mustar, -1, -2, -rhostar, 0;...
                0, 0, 0, 0, 1, 1;...
                0, 0, 0, 0, L, 0]);

            % Make PRr
            PRr=expm((log(rlayer)-log(rlayer2))*A)*PRr;
            if B>=rp
                % Make PBr
                if B<rlayer2
                    if B>rlayer
                        PBr=expm((log(rlayer)-log(B))*A)*PBr;
                    else
                        PBr=expm((log(rlayer)-log(rlayer2))*A)*PBr;
                    end
                end
            end
        end
    end

else
warning('Input variable ''precise'' should be logical 0 or 1')    
end

