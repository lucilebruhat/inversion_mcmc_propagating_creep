function [u,uc] = ve_flow_varcreep_tvec_cycle(t,tR,Nn,zco,sdot,s_rate,H,D,xobs,fault_type)
% Code to compute the surface deformation caused by viscoelasti flow
% induced by deep creep
%
% inputs: 
% - t in years, time of measurement (vector!)
% - T in years, earthquake recurrence time
% - tR in years, relaxation time in years
% - Nn, upper index fror infinite summation
% - K, number of earthquakes considered
% - D, locking depth in km
% - H, depth of visocelastic layer, in km
% - xobs, location of surface observation points in km
% - fault_type, strike slip fault orientation
% output:
% - u, surface rates, same units than v_inf
%
% Lucile Bruhat
% November 2017

% Compute Fn (which is space dependent here)
Fn = zeros(length(xobs),Nn,length(zco));

if strcmp(fault_type,'right_lat')
    
    for i = 1:length(zco)
        for n = 1:Nn % set up for right lateral strike-slip fault!!!!
            Fn(:,n,i)= -atan((2*n*H+zco(i))./xobs) - atan((zco(i)-2*n*H)./xobs);
        end
    end
    
elseif strcmp(fault_type,'left_lat')
    
    for i = 1:length(zco)
        for n = 1:Nn % set up for right lateral strike-slip fault!!!!
            Fn(:,n,i)= atan((2*n*H+zco(i))./xobs) + atan((zco(i)-2*n*H)./xobs);
        end
    end
    
end

Ds = zeros(length(zco),length(t));
for k = 1:length(t)
    for i = 1:length(zco)-1
            Ds(i,k) = sdot(i,k)-s_rate;
    end
end
    
% Obtain surface velocity due to coseismic slip delta_i between two depths
vi = zeros(length(xobs),length(t),length(zco));

for k = 1:length(t)
    tu = t(k);
    for i = 1:length(zco)-1
        % Set up for right-lateral fault
        Sumv =0;
        % loop over n
            tp = t;
            if k > 1
                for n = 1:Nn
                    pv = (Fn(:,n,i+1) - Fn(:,n,i))/factorial(n-1)*...
                        trapz(tp(1:k),((tu - tp(1:k))/tR).^(n-1).*exp(-(tu-tp(1:k))/tR).*Ds(i,1:k));
                    Sumv = Sumv + pv;
                end
            else
                for n = 1:Nn
                    pv = (Fn(:,n,i+1) - Fn(:,n,i))/factorial(n-1)*...
                        ((tu - tp(1:k))/tR).^(n-1).*exp(-(tu-tp(1:k))/tR).*Ds(i,1:k);
                    Sumv = Sumv + pv;
                end
            end
            vi(:,k,i) = Sumv/pi/tR;
    end
end

% sum to get surface deformation in m/yr
us = sum(vi,3);

%% Viscoelastic deformation due to constant creep between D and H
% Corresponds to equation 12.31

Fn = zeros(length(xobs),Nn);

if strcmp(fault_type,'right_lat')
    for n = 1:Nn % set up for right lateral  strike-slip fault!!!!
        Fn(:,n)= -atan((2*n*H+D)./xobs) - atan((D-2*n*H)./xobs);
    end  
    uc = s_rate/pi*sum(Fn,2);
    %uc = (s_rate/pi*(sum(Fn,2)-atan(xobs'./D)+atan(xobs'./H)));
elseif strcmp(fault_type,'left_lat')   
    for n = 1:Nn % set up for right lateral  strike-slip fault!!!!
        Fn(:,n)= atan((2*n*H+D)./xobs) + atan((D-2*n*H)./xobs);
    end
    uc = s_rate/pi*sum(Fn,2);
    %uc = (s_rate/pi*(sum(Fn,2)+atan(xobs'./D)-atan(xobs'./H)))
end

u = zeros(size(us));
for k = 1:length(t)
    u(:,k) = us(:,k)+uc;
end



end

