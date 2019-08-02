function u = ve_cycle_varsli(t,tR,T,Nn,K,zco,slipco,H,xobs,fault_type)
%
% inputs: 
% - t in years, time tof measurement
% - T in years, earthquake recurrence time
% - tR in years, relaxation time in years
% - Nn, upper index fror infinite summation
% - K, number of earthquakes considered
% - zco, region on the fault of coseismic slip in km
% - slipco, slip distrbution along zco in m
% - H, depth of visocelastic layer, in km
% - xobs, location of surface observation points in km
% output:
% - u, surface rates, in m/yr
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

% Compute loop over K EQ
SumEQ=zeros(1,Nn);
for n = 1:Nn
    Sumv1 = 0;
    % look over EQ numbers
    for k=0:K
        p1=exp(-k*T/tR)*(((t+k*T)/tR)^(n-1));
        Sumv1=Sumv1+p1;
    end
    SumEQ(n) = Sumv1;
end

% Obtain surface velocity due to coseismic slip delta_i between two depths
vi = zeros(length(xobs),length(zco));
for i = 1:length(zco)-1
    % Set up for right-lateral fault    
    Sumv =0;
    % loop over n
    for n = 1:Nn
        pv = SumEQ(n)*(Fn(:,n,i+1) - Fn(:,n,i))/factorial(n-1);
        Sumv = Sumv+pv;
    end
    vi(:,i) = exp(-t/tR)*Sumv/(pi*tR);  
end

% sum to get surface deformation in m/yr, convert in mm/yr
u = (vi*slipco);

end





