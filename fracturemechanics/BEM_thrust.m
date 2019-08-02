function [Gs, Gn] = BEM_thrust(x, dip, mu, nu)
% [Gs, Gn] = BEM_thrust(x, dip, mu, nu);
%
% Function to use for BEM calculation of dipping fault
%
%Input:
%  mu   shear modulus = G
%  nu   Poisson's ratio
%  x    vector of fault positions (starting at updip end and going down).
%  dip  dip of fault in degrees, dip is in positive x direction
%
%Output:
%    Gs  Shear stress at midpoint of elements due to unit slip
%           such that \Delta tau = Gs* slip;  
%    Gn  Normal stress at midpoint due to unit slip
%    Last column is effect of semi-infinite dislocation at depth
%           eg Gs = N x (N+1)
%
% Sign convention:
%           Positive slip is thrust
%           Positive shear stress is increase in driving stress
%           Positive normal stress is compression
%

%  Notes in dislocation solution, z is coordinate of dislocation, t is
%  field point

dip = dip*pi/180;

N = length(x);
deltaX = x(2)-x(1);
Gs = zeros(N-1,N);
Gn = zeros(N-1,N);

% vector of mid-points where stress is computed
t = x(1:N-1) + deltaX/2;

for i = 1:N-1
    %upper dislocation
    z = x(i);
    Ts1 = (-2*z*(t + z)*mu.*(t.^4 - 3*t.^3*z + 7*t.^2*z.^2 -  ...
        3*t*z.^3 + z.^4 + (t.^4 - 5*t.^3*z + 4*t.^2*z.^2 - 5*t*z.^3 + z.^4).*  ...
        cos(2*dip) + t.^2*z.^2*cos(4*dip))*sin(dip).^2)./  ...
       (pi*(t - z)*(-1 + nu).*(t.^2 + z.^2 - 2*t*z*cos(2*dip)).^3);
   
    Tn1 = (4*z.^2*mu*cos(dip)*(-3*t.^2*z + z.^3 + 2*t.^3*cos(2*dip)).*  ...
      sin(dip).^3)./(pi*(-1 + nu)*(t.^2 + z.^2 - 2*t*z*cos(2*dip)).^3);

    %lower dislocation
    z = x(i+1);
    Ts2 = (-2*z*(t + z)*mu.*(t.^4 - 3*t.^3*z + 7*t.^2*z.^2 -  ...
        3*t*z.^3 + z.^4 + (t.^4 - 5*t.^3*z + 4*t.^2*z.^2 - 5*t*z.^3 + z.^4).*  ...
        cos(2*dip) + t.^2*z.^2*cos(4*dip))*sin(dip).^2)./  ...
       (pi*(t - z)*(-1 + nu).*(t.^2 + z.^2 - 2*t*z*cos(2*dip)).^3);

    Tn2 = (4*z.^2*mu*cos(dip)*(-3*t.^2*z + z.^3 + 2*t.^3*cos(2*dip)).*  ...
      sin(dip).^3)./(pi*(-1 + nu)*(t.^2 + z.^2 - 2*t*z*cos(2*dip)).^3);

    Gs(:,i) = -(Ts1 - Ts2);   % minus sign makes stress increase outside slip zone
    Gn(:,i) = -(Tn1 - Tn2);
    
    
end

% Append effect of deep semi-infinite dislocation

z = x(N);
    Ts1 = (-2*z*(t + z)*mu.*(t.^4 - 3*t.^3*z + 7*t.^2*z.^2 -  ...
        3*t*z.^3 + z.^4 + (t.^4 - 5*t.^3*z + 4*t.^2*z.^2 - 5*t*z.^3 + z.^4).*  ...
        cos(2*dip) + t.^2*z.^2*cos(4*dip))*sin(dip).^2)./  ...
       (pi*(t - z)*(-1 + nu).*(t.^2 + z.^2 - 2*t*z*cos(2*dip)).^3);
   
    Tn1 = (4*z.^2*mu*cos(dip)*(-3*t.^2*z + z.^3 + 2*t.^3*cos(2*dip)).*  ...
      sin(dip).^3)./(pi*(-1 + nu)*(t.^2 + z.^2 - 2*t*z*cos(2*dip)).^3);
  
Gs(:,N) = -Ts1;

Gn(:,N) = -Tn1;
