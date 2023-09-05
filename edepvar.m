clc;
clear all;
close all;

y=0:0.5e-3:30e-3;
x=0:0.5:4;
[x1 y1]=meshgrid(x, y);
for i=1:1:size(x1,1)
    for j=1:1:size(x1,2)
        [p(i,j) ic(i,j)]=comb(y1(i,j),x1(i,j));
    end
end
