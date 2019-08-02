function Tn = virtualdislocation(t,T,tR,vinf,n)

sumv = 0;
for i = 1:20
    k = i-1;
    sum = exp(-k*T/tR)*((t+k*T)./tR).^(n-1);
    sumv = sumv+sum;
end
   
Tn = vinf*T/tR.*exp(-t/tR)/(factorial(n-1)).*sumv;

