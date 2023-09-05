clc;
clear all;
close all;
tspan=0:0.1e-9:40e-9;
options=odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-6]);
ar=500e-8;h=6.6262e-34;
dn=75;
q=1.6e-19;
tb=1e-9;
ld=sqrt(tb*dn);
dbw=250e-7;
go=3600;ns=0.26e18;ng=2.1e18;
j=1;i=1;%ie=60e-3;%vcb=-0.1;
%for vcb=0:1:2
vcb=2;ie=30e-3;
    j=1;     
%for ie=0e-3:0.5e-3:200e-3
vcb1(i,j)=vcb;
ie1(i,j)=ie;
%vcb=0;
[t, y ]=ode45(@carriersoln31,tspan,[0;0;0],options,ie,vcb,ld);
np1(i,j)=y(size(y,1),3);
np2(i,j)=y(size(y,1),2);
np3(i,j)=y(size(y,1),1);
ic(i,j)=colcurr1(np3(i,j),np1(i,j),vcb,ld);
ga(i,j)=go*log((np2(i,j)+ns)/(ng+ns));
dga(i,j)=go/(np2(i,j)+ns);
a(i,j)=ic(i,j)/ie1(i,j);
ib(i,j)=ie-ic(i,j);
p(i,j)=0.34*0.782e10*(26.19+5)*h*2.30e14*np1(i,j)*(7.5e-12/0.033);
j=j+1;
%end
i=i+1;
%end
%{

je=ie/ar;


comp1=(je/cosh(dbw/(2*ld)));

comp2=np3(j)*((sinh(dbw/(2*ld)))/(cosh(dbw/(2*ld))));

comp3=np3(j)*((cosh(dbw/(2*ld)))/(sinh(dbw/(2*ld))));

%comp4=0;
comp4=nw/sinh(dbw/(2*ld));

jvs=comp1-((q*(dn/ld))*(comp2+comp3-comp4));

ivs=jvs*ar;

jc=q*(dn/ld)*(1/sinh(dbw/(2*ld)))*(nvs - nw*cosh(dbw/(2*ld)));

ic=jc*ar;

%}

