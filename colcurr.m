function ic=colcurr(nvs,s,vcb,ld)

d=3*5e-7;%well thickness
eps=3*0.011; %optical confinement factor;

db=10e-7;%barrier thickness
w=2e-4;%stripe width
dbw=250e-7;%total base width
go=3600;%material gain
a1=5;%internal loss;
eoth=0.5e-17;%gain compression factor
beff=1.55e-10;%recombination coeff
c=1e-6;%spontaneous emission coeff'
tb=1e-9;%carrier lifetime in base
tcap=1e-12;%capture time
tesc=1e-9;%escape time
%vp=3e10/3.4850;
vg=0.782e10;
ng=2.1e18;%transparency carrier density
tp=4.1e-12;%photon lifetime
ns=0.26e18;%%fitting parameter value
q=1.6e-19;%electron charge
ar=500e-8;%area
dn=75;
na=1e18;
nd=1e18;
ni=1e10;
t=300;
k=1.38e-19;
m=9.109e-31;
me=0.048*m;
mh=0.33*m;
h=4.135e-15;
hb=h/(2*pi);
w1=2*pi*((3e10)/1.5e-4);
w0=2*pi*((3e10)/1.3e-4);
vp=2;
ip=5e-3;
mu=power(((1/me)+(1/mh)),-1);

per=8.854e-14; % electrical permittivity of free space
vbi=(k*(t/q))*(log(na*(nd/(ni^2))));
dbc=power((2*(per/q)*(vbi+vcb)*((na+nd)/(na*nd))),0.5);
f=(vbi+vcb)/dbc;
thef=power(((q^2)*(f^2)/(2*mu*hb)),1/3);
bfk=(w1-w0)/thef;

afk=((1e4*(thef^0.5))/w0)*(-bfk*power(abs(airy(0,bfk)),2)+power(abs(airy(1,bfk)),2));

nw=(0.8*ld*(d/dn)*(0.4/eps)*vg*afk*s)+((ld/(q*dn))*(ip/vp)*(vcb/ar)* exp(1+(vcb/vp)));

jc=(q*(dn/ld)*(1/sinh(dbw/(2*ld))))*(nvs-(nw*cosh(dbw/(2*ld))));

ic=jc*ar;
end