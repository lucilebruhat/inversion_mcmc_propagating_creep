function [D,Ddot,g,gp,fn,Fn,hn,tn] = ingredient_crack_fast_v3(xi,param)
% Ingredients to make the crack solution in a half-space
% Lucile Bruhat
% Modified from ingredient_crack_fast.m
% Created: August 2018
% Last modified: August 2018

% "Fast" because uses the analyticial definition of Chebyshev polynomials
% v3 because takes into account the viscoelastic relaxation with equal
% weight of T_1 and T_2 in equation 12.21 of Paul's book

N = param.N;
%if N==1; N = 2; end

v_inf = param.v_inf;% v_inf in m/s;
t = param.t;% t in s
T = param.T; % Earthquake recurrence time in s
tR = param.tR; % Relaxation time in s

m = length(xi);

% D & Ddot

% Take into account the viscoelastic relaxation
D = 1/2*virtualdislocation_disp(t,T,tR,v_inf,1) + 1/2*virtualdislocation_disp(t,T,tR,v_inf,2);
Ddot = 1/2*virtualdislocation(t,T,tR,v_inf,1) + 1/2*virtualdislocation(t,T,tR,v_inf,2);

% D = v_inf*t;
% Ddot = v_inf;

% g and gp
g = (xi.*sqrt(1-xi.^2) + asin(xi) + pi/2)/pi;
gp = 2*sqrt(1-xi.^2)/pi;

% Definition for Chebyshev polynomials 
% Do not use the matlab function for Chebyshev polynomials -> too slow
% T(n) = cos(n*acos(xi))
% U(n) = sin((n+1)*acos(xi))./sin(acos(xi))
% U(n-1) = sin(n*acos(xi))./sin(acos(xi))
% U(n-2) = sin((n-1)*acos(xi))./sin(acos(xi))

% fn
% fn(:,n)= sqrt(1-xi.^2).*(chebyshevU(n,xi)/(n+1)-chebyshevU(n-2,xi)/(n-1));
fn = zeros(m,N);
for n = 2:N
   fn(:,n)= sqrt(1-xi.^2).*(sin((n+1)*acos(xi))./sin(acos(xi))/(n+1)- ...
            sin((n-1)*acos(xi))./sin(acos(xi))/(n-1));
end
fn = fn(:,2:end);% We don't need to invert for c1 -> remove first column
fn(end,:) = zeros(1,N-1);


% Fn = fn + (1-xi)*fnp
% Fn(:,n)= sqrt(1-xi.^2).*(chebyshevU(n,xi)/(n+1)-chebyshevU(n-2,xi)/(n-1)) + ...
%            (1-xi).*(2*sqrt(1-xi.^2).*chebyshevU(n-1,xi));
Fn = zeros(m,N);
for n = 2:N
   Fn(:,n)= sqrt(1-xi.^2).*( sin((n+1)*acos(xi))./sin(acos(xi))/(n+1)-...
            sin((n-1)*acos(xi))./sin(acos(xi))/(n-1)) + ...
            (1-xi).*(2*sqrt(1-xi.^2).*sin(n*acos(xi))./sin(acos(xi)));
end
Fn = Fn(:,2:end);% We don't need to invert for c1 -> remove first column
Fn(end,:) = zeros(1,N-1);

% hn
% hn(:,n)= n*chebyshevU(n-1,xi);
hn = zeros(m,N);
for n = 2:N
   hn(:,n)= n*sin(n*acos(xi))./sin(acos(xi));
end
hn = hn(:,2:end);% We don't need to invert for c1 -> remove first column
hn(end,:) = [2:N].^2;

% tn (classic T tchebyshev polynomial chebyshevT(n,xi))
% tn(:,n)= cos(n*acos(xi));
tn = zeros(m,N);
for n = 2:N
    tn(:,n)= cos(n*acos(xi));
end
tn = tn(:,2:end);% We don't need to invert for c1: remove first column

end

