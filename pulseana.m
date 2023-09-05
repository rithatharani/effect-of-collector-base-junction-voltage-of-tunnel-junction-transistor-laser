%CE Configuration

clc;
clear all;
close all;
tspan=0:0.01e-9:40e-9;
options=odeset('RelTol',1e-4,'AbsTol',[1e-9 1e-9 1e-9]);
j=1;h=6.6262e-34;c=3e8;bef=1.55e-10;tb=1e-9;ar=500e-8;
%    %w=2*pi*10e9*t1(j);
dbw=250e-7;q=1.6e-19;
dn=75;
ld=sqrt(dn*tb);
te=dbw/(2*ld);
d=5e-7;eps=0.011;

i=1;

    j=1;
  ie=40e-3;
  vcb1=[0 7 0];
  
  tsamp=0:0.005e-9:2e-9;

  j1=1;
 for i1=1:1:length(vcb1)
     
  vcb=(vcb1(i1));
 if (j1==1)
 [t y ]=ode45(@carriersoln31,tsamp,[0;0;0],options,ie,vcb,ld);
 j1=j1+1;
 else
 [t y ]=ode45(@carriersoln31,tsamp,[np3(i-1,size(np3(i-1,:),2));np2(i-1,size(np2(i-1,:),2));np1(i-1,size(np1(i-1,:),2))],options,ie,vcb,ld);
 end
np1(i,:)=y(:,3); % photon density 
np2(i,:)=y(:,2); % QW electron density
np3(i,:)=y(:,1); % VS electron density 


for j=1:1:length(np1(i,:))
% Virtual State Current Density
i1(i,j)=ie;
ic(i,j)=colcurr1(np3(i,j),np1(i,j),vcb,ld);
ib(i,j)=i1(i,j)-ic(i,j);

%bet(i,j)=ic(i,j)/i1(j);

%optical Power
p(i,j)=0.34*0.782e10*(26.19+5)*h*2.30e14*np1(i,j)*(7.5e-12/0.033);
end
i=i+1;
 end
 p1=reshape(p.',1,[]);
 np21=reshape(np2.',1,[]);
   np11=reshape(np1.',1,[]);
   np31=reshape(np3.',1,[]);
  
 t1=0:0.005e-9:0.005e-9*(length(p1)-1);

 i11=zeros(1,length(t1));len=0;
 for i2=1:1:length(vcb1)
 i11(1+len:len+(length(t1)/length(vcb1)))=vcb1(i2);
 len=len+(length(t1)/length(vcb1));
 end
 
 figure(1);
 subplot(2,1,1);
  plot(t1(length(p(1,:)):length(t1)),i11(length(p(1,:)):length(p1)));
  subplot(2,1,2);
    plot(t1(length(p(1,:)):length(t1)),p1(length(p(1,:)):length(p1)));
    
    
 figure(2);
   plot(t1,i11);
   figure (3);
      plot(t1,p1);
     figure (4);
     subplot(3,1,2)
  plot(t1,np21);
  subplot(3,1,1);
  plot(t1,np31);
subplot(3,1,3)
    plot(t1,np11);
