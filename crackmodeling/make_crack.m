function [dtau,dtau_rate,slip,slip_rate] = make_crack(xi,cn,cndot,adot,param)
% Crack solution in a half-space
% Need to compute the "ingredients" fist from ingredient_crack_fast
% Lucile Bruhat
% Created: October 2016
% Last modified: July 19 (corrected factor 2 in stress definition)
% June 22, 2017 (changed v_inf*t to D to be more general)

% All parameters have to be in common units (Pa,m,s)

a = param.a;
mu = param.mu;
nu = param.nu;
D = param.D;
Ddot = param.Ddot;
g = param.g;
gp = param.gp;
fn = param.fn;
Fn = param.Fn;
hn = param.hn;
tn = param.tn;
N = param.N;

% antiplane vs in-plane deformation
if (param.mode == 2) 
    mu_s = mu/(1-nu);
elseif (param.mode == 3)
    mu_s = mu;
else 
    warning('Unknown fracture mode')  
end

if N > 1
    % reconstruct slip rate
    slip_rate = Ddot*g + a/2*fn*cndot +adot*((1-xi).*D.*gp/a+Fn*cn/2);

    % reconstruct slip at time t
    slip = D*g + a/2*fn*cn;
    
    % reconstruct stress drop
    dtau = mu_s*(2*D/a/pi*xi + tn*cn);
    
    % reconstruct stress drop change in Pa/s
    dtau_rate = mu_s*(2*Ddot/a/pi*xi + tn*cndot + adot/a*((1-2*xi)*2*D/a/pi + (1-xi).*(hn*cn)));

else   % If N = 1, all cn and cndot are zeros
    
    % reconstruct slip rate
    slip_rate = Ddot*g  +adot*(1-xi).*D.*gp/a;
    
    % reconstruct slip at time t
    slip = D*g;
    
    % reconstruct stress drop
    dtau = mu_s*2*D/a/pi*xi;
    
    % reconstruct stress drop change in Pa/s
    dtau_rate = mu_s*(2*Ddot/a/pi*xi + adot/a*(1-2*xi)*2*D/a/pi);
    
    %acceleration
    gpp = -2*xi./(pi*sqrt(1-xi.^2));
    slip_acc = adot*2*(1-xi).*gp*Ddot/a + adot^2/(a^2)*D*(1-xi).*((1-xi).*gpp - 2.*gp) ; % to check
    
end