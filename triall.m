clc;
clear all;
close all;

nvs=1e18;
wb=250e-7;
dn=75;
tb=1e-9;
ld=sqrt(dn*tb);
ie=1e-4;
ar=500e-8;
je=ie/ar;
q=1.6e-19;

den1=2*sinh((wb/(2*ld)));

comp11=(1/den1)*(nvs*exp((wb/(2*ld))))-((ld*je)/(q*dn))
comp12=(1/den1)*((ld*je)/(q*dn))-(nvs*exp(((-wb)/(2*ld))));

den2=2*cosh((wb/(2*ld)));
comp21=(1/den2)*(nvs*exp((wb/(2*ld))))-((ld*je)/(q*dn));
comp22=(1/den2)*(nvs*exp(((-wb)/(2*ld)))+((ld*je)/(q*dn)));

i=1;
for x=-0.1:0.001:0
n1(i)=comp11*exp(x/ld)+comp12*exp((-x)/ld);

n2(i)=comp21*exp(x/ld)+comp22*exp((-x)/ld);
x1(i)=x;
i=i+1;
end


