
%CB Configuration
clc;
clear all;
close all;

options=odeset('RelTol',1e-4,'AbsTol',[1e-9 1e-9 1e-9]);
h=6.6262e-34;bef=1.55e-10;ar=500e-8;

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
ieb1=34e-3*[1.5];
for j4=1:1:length(ieb1)
ieb=ieb1(j4);
vcb=2.0;
carr=freqsub(ieb,vcb);
pb=carr(4);
j2=1;j3=1;
for f=1e8:1e8:10e9
%f=2.35e9;%input Frequency
%f1=2.25e9;

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
    tspan=0:0.1e-9:9e-9;
end

for i4=1:length(t1)
  i2(j3,i4)=ieb +(1.5e-3)*sin(2*pi*f*t1(i4));%+(1.5e-3)*sin(2*pi*f1*t1(i4));%Input Emitter Current
end
j=1;j1=1;
for i3=1:length(i2(j3,:))
if j==1
    %[ti yi ]=ode45(@carriersoln31,tspan,[0;0;0],options,i2(j3,i3),vcb,ld);
[t y ]=ode45(@carriersoln31,tspan,[carr(1);carr(2);carr(3)],options,i2(j3,i3),vcb,ld);
else
    [t y ]=ode45(@carriersoln31,tspan,[np3(j3,j1);np2(j3,j1);np1(j3,j1)],options,i2(j3,i3),vcb,ld);
    j1=j1+1;
end
    
np1(j3,j)=y(length(y),3);
np2(j3,j)=y(length(y),2);
np3(j3,j)=y(length(y),1);
i1(j3,j)=i2(i3);
ic(j3,j)=colcurr1(np3(j3,j),np1(j3,j),vcb,ld);
ib(j3,j)=i1(j3,j)-ic(j3,j);
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
db(j4,j2)=10*log10(po(j2)/(pb));
fre(j4,j2)=f;
j2=j2+1;
j3=j3+1;
end
end


%{
j2=j2-1;
j3=j3-1;

plot(t1,p(j3,:));
title('optical power output');
xlabel('time');
ylabel('power in watts');
figure;
plot(t1,np1(j3,:));
title('optical output');
xlabel('time');
ylabel('photon density/cm3 ');
figure;
plot(t1,np2(j3,:));
title('electron density');
xlabel('time');
ylabel('electron density/cm3 ');
figure;
plot(t1,i2(j3,:));
title('Input emitter signal');
xlabel('time');
ylabel('Amplitude in amperes ');
figure;
plot(t1,ic(j3,:));
title(' collector signal');
xlabel('time');
ylabel('Amplitude in amperes ');
figure;
plot(t1,ib(j3,:));
title(' base signal');
xlabel('time');
ylabel('Amplitude in amperes ');
%}

ma=figure;
semilogx(fre,db);
title('Frequency plot');
xlabel('frequency');
ylabel('magnitude in dB ');
%saveas(ma,'F:\project\optics\review\Phase II\paper\CB\3qw\425ma.fig')
%plot(t1,p);

%figure;

%plot(t1,i1*1e3);



 
 
 