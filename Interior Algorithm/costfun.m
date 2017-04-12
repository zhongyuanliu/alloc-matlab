function y=costfun(x,x0)
w1=100;
w4=50;
w2=10;
w5=5;
w3=1;
w6=.5;
% syms w1 w2 w3 w4 w5 w6 real
y=w1*(x(1)-x0(1))^2+w2*(x(2)-x0(2))^2+w3*(x(3)-x0(3))^2 ...
+w4*(x(4)-x0(4))^2+w5*(x(5)-x0(5))^2+w6*(x(6)-x0(6))^2;

end