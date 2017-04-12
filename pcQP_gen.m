function [ x_stop, y_stop, z_stop, s_stop, k ] = pcQP_gen (x, y, z, s,G, g ,C, d ,A, b)
%----------------------------------------------------------
% pcQP gen .m
%
% This f u n c t i on s o l v e s a QP problem o f the form
% min 0 . 5 x 'Gx + g ' x
% s . t . A' x = b
% C' x >= d
% us ing the Pr edi c tor-Cor r e c tor (PC) method .
%
% Input :
% x : s t a r t i n g point f o r the ve c tor x (n x 1 ve c tor )
% y : s t a r t i n g point f o r the ve c tor x (nA x 1 ve c tor )
% z : s t a r t i n g point f o r the ve c tor x (nC x 1 ve c tor )
% s : s t a r t i n g point f o r the s lack-ve c tor x (nC x 1 ve c tor )
% G: He s s ian (n x n matrix )
% g : (n x 1 ve c tor )
% C: l e f t hand s i d e o f i n e q u a l i t y c o n s t r a i n t s (n x nC matrix )
% d : r i gh t hand s i d e o f i n e q u a l i t y c o n s t r a i n t s (nC x 1 ve c tor )
% A: l e f t hand s i d e o f e qu a l i t y c o n s t r a i n t s (n x nA matrix )
% b : r i gh t hand s i d e o f e qu a l i t y c o n s t r a i n t s (nA x 1 ve c tor )
% where nC and nA ar e the numbers o f i n e q u a l i t y and e qu a l i t y
% c o n s t r a i n t s .
% Output :
% x_stop : s o l u t i on x
% y_stop
% z s t op
% s s t op
% k : number o f i t e r a t i o n s used
%
% Thomas Reslow Kr¡§uth , s021898
%----------------------------------------------------------
%dampening f a c t o r eta
eta = 0.95 ;
%r e s i d u a l s ar e computed
[mA,nA] = size (A) ;
[mC, nC] = size (C) ;
e = ones (nC, 1) ;
rL = G*x + g - A*y - C*z ;
rA = -A' * x + b ;
rC = -C' * x + s + d ;
rsz= s .* z ;
mu = sum( z .* s ) /nC;
%k : number o f i t e r a t i o n s , e p s i l on : t o l e r a n c e s
k = 0 ;
maxk = 200;
eps_L = 1e-10; eps_A = 1e-10; eps_C = 1e-10; eps_mu = 1e-10;
while ( k<=maxk && norm( rL )>=eps_L && norm( rA)>=eps_A && norm( rC)>=eps_C && abs (mu)>=eps_mu )
%Solve system with a Newton-l i k e method/ Fa c t o r i z a r i on
lhs = [G,-A,-C;-A' , spalloc(nA, nA, 0 ) , spalloc(nA, nC, 0 ) ;-C' , spalloc(nC,nA, 0 ) , sparse(-diag ( s ./ z ) ) ] ;
[L,D,P] = ldl ( lhs ) ;
rhs = [-rL;-rA;-rC+rsz./ z ] ;
dxyz_a = P*(L' \ (D\(L\(P' * rhs ) ) ) ) ;
dx_a = dxyz_a ( 1 : length ( x ) ) ;
dy_a = dxyz_a ( length ( x )+1: length ( x )+length ( y ) ) ;
dz_a = dxyz_a ( length ( x )+length ( y )+1: length ( x )+length ( y )+length ( z ) ) ;
ds_a = -(( rsz+s .* dz_a ) ./ z ) ;
%Compute a l p h a a f f
alpha_a = 1 ;
idx_z = find( dz_a <0) ;
if ( isempty ( idx_z )==0)
alpha_a = min ( alpha_a ,min(-z( idx_z ) ./ dz_a( idx_z ) ) ) ;
end
idx_s = find( ds_a <0) ;
if ( isempty ( idx_s )==0)
alpha_a = min ( alpha_a ,min(-s ( idx_s ) ./ ds_a ( idx_s ) ) ) ;
end
%Compute the a f f i n e d u a l i t y gap
mu_a = ( ( z+alpha_a * dz_a )' * ( s+alpha_a * ds_a ) ) /nC;
%Compute the c e n t e r i n g parameter
sigma = (mu_a/mu)^3 ;
%Solve system
rsz= rsz+ ds_a .* dz_a - sigma * mu * e ;
rhs = [-rL;-rA;-rC+rsz./ z ] ;
dxyz = P*(L' \ (D\(L\(P' * rhs ) ) ) ) ;
dx = dxyz ( 1 : length ( x ) ) ;
dy = dxyz ( length ( x )+1: length ( x )+length ( y ) ) ;
dz = dxyz ( length ( x )+length ( y )+1: length ( x )+length ( y )+length ( z ) ) ;
ds = -(( rsz + s .* dz ) ./ z ) ;
%Compute alpha
alpha = 1 ;
idx_z = find( dz<0) ;
if ( isempty ( idx_z )==0)
alpha = min( alpha ,min(-z( idx_z ) ./ dz( idx_z ) ) ) ;
end
idx_s = find( ds<0) ;
if ( isempty ( idx_s )==0)
alpha = min ( alpha ,min(-s( idx_s )./ ds( idx_s ) ) ) ;
end
%Update x , z , s
x = x + eta * alpha*dx ;
y = y + eta * alpha*dy ;
z = z + eta * alpha*dz ;
s = s + eta * alpha*ds ;
k = k+1;
%Update rhs
rL = G*x + g - A*y - C*z ;
rA = -A' * x + b ;
rC = -C' * x + s + d ;
rsz= s .* z ;
mu = sum( z .* s ) /nC;
end
%Output
x_stop = x ;
y_stop = y ;
z_stop = z ;
s_stop = s ;