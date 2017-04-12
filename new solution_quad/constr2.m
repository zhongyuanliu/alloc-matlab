function [y]=constr2(x,n_2)
%n_2=n4-n1
% syms n_2 real
n2=sin(x(1))*x(4);
n3=n2+sin(x(1)-(pi-x(2)))*x(5);
n4_=n3+ sin(x(1)-(pi-x(2))-(pi-x(3)))*x(6);
% c=[-((Rbase_n*[n3(1)-n4_(1) 0 n3(2)-n4_(2)]')'*[0 0 1]')
%     -((Rbase_n*[n3(1)-n4_(1) 0 n3(2)-n4_(2)]')'*V5n)];%slow
% c=[];
% c=[x(1)-(pi-x(2))-(pi-x(3))-10
%    -80/180*pi-(x(1)-(pi-x(2))-(pi-x(3)))];
y=n4_-n_2;