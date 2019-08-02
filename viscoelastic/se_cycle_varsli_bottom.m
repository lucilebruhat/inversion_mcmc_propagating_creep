function s = se_cycle_varsli_bottom(t,tR,T,mu,Nn,K,zslip,slip,H,zobs,fault_type)
%
% inputs:
% - t in years, time tof measurement
% - T in years, earthquake recurrence time
% - tR in years, relaxation time in years
% - mu, shear modulus in Pa
% - Nn, upper index fror infinite summation
% - K, number of earthquakes considered
% - zco, region on the fault of coseismic slip in km
% - slip, slip distrbution along zco in m
% - H, depth of visocelastic layer, in km
% - zobs, depth of observation points in km (have to be >H)
% output:
% - s, stress on the fault
%
% Lucile Bruhat
% June 2018

zslip(1) = [];slip(1) = [];
zslip(end) = [];slip(end) = [];
zobs(1) = [];

% zobs has to be larger than H
if min(zobs)<H
    warning('Observation points have to be in the viscoelastic medium')
    return
end

% Compute Wn (which is space dependent here)
Wn = zeros(Nn,length(zobs),length(zslip));

if strcmp(fault_type,'right_lat')
    
    for i = 1:length(zobs)
        for j = 1:length(zslip)
            for n = 1:Nn % set up for right lateral strike-slip fault!!!!
                Wn(n,i,j)= 1/(zobs(i)+2*n*H+zslip(j)) - 1/(zobs(i)+2*n*H-zslip(j));
            end
        end
    end
     
elseif strcmp(fault_type,'left_lat')
    
    for i = 1:length(zobs)
        for j = 1:length(zslip)
            for n = 1:Nn
                Wn(n,i,j)= - 1/(zobs(i)+2*n*H+zslip(j)) + 1/(zobs(i)+2*n*H-zslip(j));            
            end
        end
    end
    
end

% Compute loop over K EQ
SumEQ=zeros(1,Nn);
for n = 1:Nn
    Sumv1 = 0;
    % look over EQ numbers
    for k=0:K
        p1=exp(-k*T/tR)*((((t+k*T)/tR)^(n-1))-1);
        Sumv1=Sumv1+p1;
    end
    SumEQ(n) = Sumv1;
end



% Obtain surface velocity due to coseismic slip delta_i between two depths
si = zeros(length(zobs),length(zslip));

for j = 1:length(zslip)-1
    % Set up for right-lateral fault
    Sums =0;
    % loop over n
    for n = 1:Nn
        ps = SumEQ(n)*(Wn(n,:,j+1) - Wn(n,:,j))/factorial(n-1);
        Sums = Sums+ps;
    end
    
    si(:,j) = mu*exp(-t/tR)*Sums/(pi*tR);
end 

s = (si*slip);

end





