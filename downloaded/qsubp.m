%求解子问题
function [x,lambda]=qsubp(H,c,Ae,be)
ginvH=pinv(H);
[m,n]=size(Ae);
if m>0
    rb=Ae*ginvH*c+be;
    lambda=pinv(Ae*ginvH*Ae')*rb;
    x=ginvH*(Ae'*lambda-c);
else
    x=-ginvH*c;
    lambda=0;
end