function y=Jh(x)
y=[ x(5)*sin(x(1) + x(2)) - x(6)*sin(x(1) + x(2) + x(3)) - x(4)*sin(x(1)), x(5)*sin(x(1) + x(2)) - x(6)*sin(x(1) + x(2) + x(3)), -x(6)*sin(x(1) + x(2) + x(3)), cos(x(1)), -cos(x(1) + x(2)), cos(x(1) + x(2) + x(3))
 x(6)*cos(x(1) + x(2) + x(3)) - x(5)*cos(x(1) + x(2)) + x(4)*cos(x(1)), x(6)*cos(x(1) + x(2) + x(3)) - x(5)*cos(x(1) + x(2)),  x(6)*cos(x(1) + x(2) + x(3)), sin(x(1)), -sin(x(1) + x(2)), sin(x(1) + x(2) + x(3))];
end