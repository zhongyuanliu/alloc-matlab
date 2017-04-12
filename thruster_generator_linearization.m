function [u0,x0,A,B] = thruster_generator_linearization(num,a,b,c,N,P)
% generate linearized constraint for thrusters connector to a same
% generator
% R^2: max power of the generator
% N: circle is divided into N parts
% u: x, y thrust of each thruster, colum vector
% sum((a .* u + bi).^2 + c) <= P, P is max power
% a,b,c colum vector, a,b,c is the 2 order polynormial fit coefficient
% x == a .* u + b
% R^2 == P - sum(c)
% sum(x.^ 2) <= R^2
% output meets A*u <= B
% [Rangle_T, ~] = polyfit(rudder_table(:,1), rudder_table(:,5),N_poly);%rudder angle vs T
% polyval(Fangle_Rangle,abs(alloc_out(i).phi))
a = abs(a);
delta = 2 * pi / N;
R = sqrt(abs(P - sum(c)));%!!must be real
r = R * cos(delta / 2);
% num = numel(u);% number of variables
x0 = ones(N ^ (num - 1),num);
% num - n_loop == 1
x0 = r * recursive_loop(x0,N,delta,num - 1);

A = x0 .* repmat(a',size(x0,1),1);
B = sum(x0.^2,2) - x0 * b;
u0 = (x0 - repmat(b',size(x0,1),1)) ./ repmat(a',size(x0,1),1);
end

function y = recursive_loop(x0,N,delta,n_loop)
if n_loop == 1
    for i = 1 : N
        x0(i,1) = cos((i - 1) * delta);
        x0(i,2) = sin((i - 1) * delta);
    end
elseif n_loop > 1
    y = recursive_loop(x0,N,delta,n_loop - 1);
    for j = 1 : N
        for n = 1 : n_loop + 1
            for i = 1 : N^(n_loop - 1)
                if n < n_loop + 1
                    x0(i + N^(n_loop - 1) * (j - 1),n) = sin((j - 1) * delta) * y(i,n);
                else
                    x0(i + N^(n_loop - 1) * (j - 1),n) = cos((j - 1) * delta) * y(i,n);
                end
            end
        end
    end
end
y = x0;
end
% a = [0.137477270848675;0.200000000000000;0.100000000000000];
% b = [9.44519768238120;12;23];
% c = [-103.511759259259;-136;-80];
% P = 1000;
% N = 200;
% [u0,x0,A,B] = thruster_generator_linearization(num,a,b,c,N,P)
% plot3(u0(:,1),u0(:,2),u0(:,3))
% axis equal
% hold on
% plot3(x0(:,1),x0(:,2),x0(:,3))