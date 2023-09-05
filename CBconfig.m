%CB Configuration
clc;
clear all;
close all;
tspan=0:0.1e-9:40e-9;
options=odeset('RelTol',1e-4,'AbsTol',[1e-9 1e-9 1e-9]);
h=6.6262e-34;
dn=75;tb=1e-9;dbw=250e-7;
j3=1;j=1;

ie=29.1e-3;
vcb=0;

ld=sqrt(dn*tb);%Diffusion Length
%ld=sqrt(dn*tb);
te=dbw/(2*ld);


[t y ]=ode45(@carriersoln31,tspan,[0;0;0],options,ie,vcb,ld);
    
np1(j3,j)=y(size(y,1),3);
np2(j3,j)=y(size(y,1),2);
np3(j3,j)=y(size(y,1),1);
i1(j3,j)=ie;
ic(j3,j)=colcurr1(np3(j3,j),np1(j3,j),vcb,ld);
ib(j3,j)=i1(j3,j)-ic(j3,j);
%Optici3al Power 
p(j3,j)=0.34*0.782e10*(26.19+5)*h*2.30e14*np1(j3,j)*(7.5e-12/0.033);
plot(t,y(:,3));

