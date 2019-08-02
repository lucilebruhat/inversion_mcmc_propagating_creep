function Dn = virtualdislocation_disp(t,T,tR,vinf,n)

x = linspace(0,t,500);

Dn = trapz(x,virtualdislocation(x,T,tR,vinf,n));