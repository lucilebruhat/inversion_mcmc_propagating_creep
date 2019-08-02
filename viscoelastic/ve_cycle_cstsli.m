function u = ve_cycle_cstsli(t,tR,T,Nn,K,v_inf,H,D,xobs,fault_type)
%
% Set up for right-lateral fault
% inputs: 
% - t in years, time tof measurement
% - T in years, earthquake recurrence time
% - tR in years, relaxation time in years
% - Nn, upper index fror infinite summation
% - K, number of earthquakes considered
% - v_inf, long-term rate, usually in mm/yr
% - D, locking depth in km
% - H, depth of visocelastic layer, in km
% - xobs, location of surface observation points in km
% - fault_type, strike slip fault orientation
% output:
% - u, surface rates, same units than v_inf
%
% Lucile Bruhat
% November 2017

% Compute Fn
Fn = zeros(length(xobs),Nn);

if strcmp(fault_type,'right_lat')
    for n = 1:Nn % set up for right lateral  strike-slip fault!!!!
        Fn(:,n)= -atan((2*n*H+D)./xobs) - atan((D-2*n*H)./xobs);
    end    
elseif strcmp(fault_type,'left_lat')   
    for n = 1:Nn % set up for right lateral  strike-slip fault!!!!
        Fn(:,n)= atan((2*n*H+D)./xobs) + atan((D-2*n*H)./xobs);
    end
    
end

Sumv =0;
% loop over n
for n = 1:Nn
    Sumv1=0;
    
    % look over EQ numbers
    for k=0:K
        p1=exp(-k*T/tR)*(((t+k*T)/tR)^(n-1));
        Sumv1=Sumv1+p1;
    end
    
    pv = Sumv1*Fn(:,n)/factorial(n-1);
    Sumv = Sumv+pv;
end

u = v_inf*T/(pi*tR)*exp(-t/tR)*Sumv;

end

