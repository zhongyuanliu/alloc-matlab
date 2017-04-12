%{
H=[1 -1; -1 2];
c=[-2 ;-6];
A=[1 1; -1 2];
b=[2;2];
Aeq=[];beq=[];
VLB=[0;0];
VUB=[];
[x,z]=quadprog(H,c,A,b,Aeq,beq,VLB,VUB)
%}
%{
clc;
clear;

syms x1 x2 t s;
f=x1^2+x2^2+6*x1+2;
[x,mf]=minRosen(f,[-2 -1;-3 -5],[-4;-8],[0 0],[t s])
%}
%{
H=[2 -2 0;-2 4 0;0 0 2];
c=[0 0 1]';
A=[1 1 1;2 -1 1];
b=[4 2]';
[x,lam,fval]=qlag(H,A,b,c)
%}
clc;
clear;
%function callqpact
H=[2 -1;-1 4];
c=[-1 -10]';
Ae=[];be=[];
Ai=[-3 -2;1 0;0 1];
bi=[-6 0 0]';
x0=[0 0]';
[x,lambda,exitflag,output]=qpact(H,c,Ae,be,Ai,bi,x0)