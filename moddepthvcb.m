%CB Configuration
clc;
clear all;
close all;

options=odeset('RelTol',1e-4,'AbsTol',[1e-9 1e-9 1e-9]);
h=6.6262e-34;bef=1.55e-10;j2=1;j3=1;ar=500e-8;
%pb=4.4344e-4;
%sp=1e5;
i=1;f1=1e8;a=1;tb=1e-9;
dbw=250e-7;q=1.6e-19;
dn=75;

%{
    for i=1:54
    f1(i+1)=f1(i)+sp;
     if (f1(i+1)==10*sp)
    sp=sp*10;
     end
     end
%}
sc=0;ll=1;



%for f=1e8:1.5e8:5e9
%for vcb=1:0.5:2.5
vcb=1;
f=2.4e9;%input Frequency
%f1=2.25e9;
ie=32e-3;
w=2*pi*f;
ww=sqrt(w^2*tb^2+1);
ld=sqrt(dn*(tb/ww));%Diffusion Length
%ld=sqrt(dn*tb);
te=dbw/(2*ld);

fs=50*f;

sc=(20e-9*f);
t1=0:1/fs:(5+sc)/f;
if (f>=25e6)
tspan=0:1/(10*fs):1/fs;
else 
    tspan=0:0.1e-9:10e-9;
end

for i4=1:length(t1)
  vcb1(j3,i4)=vcb +(1)*sin(2*pi*f*t1(i4));%+(1.5e-3)*sin(2*pi*f1*t1(i4));%Input Emitter Current
end
j=1;j1=1;
for i3=1:length(vcb1(j3,:))
if j==1
    [ti yi ]=ode45(@carriersoln31,tspan,[0;0;0],options,ie,vcb1(j3,i3),ld);
[t y ]=ode45(@carriersoln31,tspan,[yi(length(yi),1);yi(length(yi),2);yi(length(yi),3)],options,ie,vcb1(j3,i3),ld);
else
    [t y ]=ode45(@carriersoln31,tspan,[np3(j3,j1);np2(j3,j1);np1(j3,j1)],options,ie,vcb1(j3,i3),ld);
    j1=j1+1;
end
    
np1(j3,j)=y(length(y),3);
np2(j3,j)=y(length(y),2);
np3(j3,j)=y(length(y),1);
vcb2(j3,j)=vcb1(i3);
ic(j3,j)=colcurr1(np3(j3,j),np1(j3,j),vcb1(j3,i3),ld);
ib(j3,j)=ie-ic(j3,j);
%Optici3al Power 
p(j3,j)=0.34*0.782e10*(26.19+5)*h*2.30e14*np1(j3,j)*(7.5e-12/0.033);
j=j+1;
end
%{
if f>=1e8
    temp1=ceil(f/1e9);
    ll=50*temp1;
    if f>=8e9
        ll=ll*8;
    end
end
%}
temp3=length(t1)-length(0:1/fs:((5+sc)-2)/f);
temp2=1;
temp4=length(t1)-temp3;
for q=temp4:length(p(j3,:))
    p1(j3,temp2)=p(j3,q);
    temp2=temp2+1;
end
po(j2)=max(p1(j3,:))-min(p1(j3,:));
%poi=power((max(i1(1,:))),2)-power((min(i1(1,:))),2);
%poi=max(poi1(1,:))-min(poi1(1,:));
%db(j2)=10*log10(po(j2)/poi);
%db(j2)=10*log10(po(j2)/(pb));
db(j2)=10*log10(po(j2)/max(p1(j3,:)));
fre(j2)=f;
j2=j2+1;
j3=j3+1;
%end

j2=j2-1;
j3=j3-1;