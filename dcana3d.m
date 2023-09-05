%CB Configuration
clc;
clear all;
close all;
tspan=0:0.1e-9:40e-9;
options=odeset('RelTol',1e-4,'AbsTol',[1e-9 1e-9 1e-9]);
h=6.6262e-34;
dn=75;tb=1e-9;dbw=250e-7;
j3=1;j=1;g0=3600;ns=0.26e18;

ie2=0e-3:0.5e-3:50e-3;
vcb2=0:0.5:3;
[ie1 vcb1]=meshgrid(ie2,vcb2);

for i=1:1:size(ie1,1)
for i1=1:1:size(ie1,2)
%vcb=0;

ld=sqrt(dn*tb);%Diffusion Length
%ld=sqrt(dn*tb);
te=dbw/(2*ld);


[t y ]=ode45(@carriersoln31,tspan,[0;0;0],options,ie1(i,i1),vcb1(i,i1),ld);
    
np1(i,i1)=y(size(y,1),3);
np2(i,i1)=y(size(y,1),2);
np3(i,i1)=y(size(y,1),1);
ie3(i,i1)=ie1(i,i1);
ic(i,i1)=colcurr1(np3(i,i1),np1(i,i1),vcb1(i,i1),ld);
ib(i,i1)=ie3(i,i1)-ic(i,i1);


%Optici3al Power 
p(i,i1)=0.34*0.782e10*(26.19+5)*h*2.30e14*np1(i,i1)*(7.5e-12/0.033);

end
    
end

for i=1:1:size(p,1)
    plot3(ie1(i,:),vcb1(i,:),p(i,:));
    hold on
end
hold off

