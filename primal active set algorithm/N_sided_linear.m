function [A,b,N] = N_sided_linear(error)
%%
N = floor(pi / (acos(1 - error))) + 1;
%     r = - R * cos(pi / N);%change <= to >=
phi = zeros(N,1);
A = zeros(N,2);
b =zeros(N,1);
for i = 0 : (N - 1)
    phi(i + 1) = (2 * i + 1) * pi / N;
    A(i + 1,1) = - cos(phi(i + 1));%change <= to >=
    A(i + 1,2) = - sin(phi(i + 1));%change <= to >=
    b(i + 1) = -(5 + i);
end
A = -A;
b = -b;
end