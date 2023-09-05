%CB Configuration
clc;
clear all;
close all;
tspan=0:0.1e-9:40e-9;
options=odeset('RelTol',1e-4,'AbsTol',[1e-9 1e-9 1e-9]);
h=6.6262e-34;
dn=75;tb=1e-9;dbw=250e-7;
j3=1;j=1;g0=3600;ns=0.26e18;ng=2.1e18;
for vcb=0:0.5:0 %ie=24e-3:1e-3:26e-3 
    j=1;
for ie=0e-3:10e-3:200e-3 % vcb=0:0.2:3
%vcb=1;

ld=sqrt(dn*tb);%Diffusion Length
%ld=sqrt(dn*tb);
te=dbw/(2*ld);


[t y ]=ode45(@carriersoln31,tspan,[0;0;0],options,ie,vcb,ld);
    
np1(j3,j)=y(size(y,1),3);
np2(j3,j)=y(size(y,1),2);
np3(j3,j)=y(size(y,1),1);
i1(j3,j)=ie;
vcb1(j3,j)=vcb;
ic(j3,j)=colcurr1(np3(j3,j),np1(j3,j),vcb,ld);
bte(j3,j)=ic(j3,j)/i1(j3,j);
ib(j3,j)=i1(j3,j)-ic(j3,j);
num=np2(j3,j)+ns;
den=ng+ns;
g(j3,j)=g0*log(num/den);
dg(j3,j)=g0/num;

%Optici3al Power 
p(j3,j)=0.34*0.782e10*(26.19+5)*h*2.30e14*np1(j3,j)*(7.5e-12/0.033);
j=j+1;
end
j3=j3+1;
end
