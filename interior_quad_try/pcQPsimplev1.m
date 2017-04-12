function [ xstop , zstop , sstop , k ] = pcQPsimplev1( x , z , s ,G, g ,A, b)
%??????????????????????????????????????????????????????????
% pcQP simplev1 .m
%
% This f u n c t i on s o l v e s a QP problem o f the form
%
% min 0 . 5 x ��Gx + g �� x
% s . t . A�� x >= b
%
% us ing the Pr edi c tor?Cor r e c tor (PC) method .
%
% Input :
% x : s t a r t i n g point f o r the ve c tor x (n x 1 ve c tor )
% z : s t a r t i n g point f o r the ve c tor x (nA x 1 ve c tor )
% s : s t a r t i n g point f o r the s lack?ve c tor x (nA x 1 ve c tor )
% G: He s s ian (n x n matrix )
% g : (n x 1 ve c tor )
% A: l e f t hand s i d e o f i n e q u a l i t y c o n s t r a i n t s (n x nA matrix )
% b : r i gh t hand s i d e o f i n e q u a l i t y c o n s t r a i n t s (nA x 1 ve c tor )
% where nC and nA ar e the numbers o f i n e q u a l i t y and e qu a l i t y
% c o n s t r a i n t s .
% Output :
% xstop : s o l u t i on x
% z s t op
% s s t op
% k : number o f i t e r a t i o n s used
%
% Thomas Reslow Kr��uth , s021898
%??????????????????????????????????????????????????????????
eta = 0.95;
%r e s i d u a l s ar e computed
[mA,nA] = size (A) ;
e = ones (nA, 1 ) ;
rL = G*x + g - A*z ;
rs = s - A' * x + b ;
rsz = s .* z ;
mu = sum( z .* s )/nA;
%k : number o f i t e r a t i o n s , e p s i l o n : t o l e r a n c e s
k = 0 ;
maxk = 200;
eps_L = 1e-16; eps_s = 1e-16; eps_mu = 1e-16;
while (k<=maxk && norm( rL )>=eps_L && norm( rs )>=eps_s && abs (mu)>=eps_mu )
%Solve system with a Newton?l i k e method/ Fa c t o r i z a r i on
G_bar = G + A*( diag ( z ./ s ) ) * A' ;
r_bar = A * ( ( rsz-z .* rs ) ./ s ) ;
g_bar = - (rL + r_bar ) ;
L = chol(G_bar , 'lower' ) ;
dx_a= L' \ (L\ g_bar ) ;
ds_a = -rs + A' * dx_a;
dz_a = -( rsz+z .* ds_a ) ./ s ;
%Compute a l p h a a f f
alpha_a = 1 ;
idx_z = find ( dz_a <0) ;
if ( isempty ( idx_z )==0)
alpha_a = min( alpha_a ,min(-z( idx_z ) ./ dz_a( idx_z ) ) ) ;
end
idx_s = find ( ds_a <0) ;
if ( isempty ( idx_s )==0)
alpha_a = min( alpha_a ,min(-s( idx-s ) ./ ds_a( idx_s ) ) ) ;
end
%Compute the a f f i n e d u a l i t y gap
mu_a = ( ( z+alpha_a * dz_a )' * ( s+alpha_a * ds_a ) ) /nA;
%Compute the c e n t e r i n g parameter
sigma = (mu_a/mu)^3 ;
%Solve system
rsz = rsz + ds_a .* dz_a - sigma * mu * e ;
r_bar = A * ( ( rsz - z .* rs ) ./ s ) ;
g_bar = -(rL + r_bar ) ;
dx = L' \ (L\ g_bar ) ;
ds = -rs + A' * dx ;
dz = -( rsz+z .* ds ) ./ s ;
%Compute alpha
alpha = 1 ;
idx_z = find(dz<0) ;
if ( isempty ( idx_z )==0)
alpha = min ( alpha ,min(-z( idx_z ) ./ dz( idx_z ) ) ) ;
end
idx_s = find ( ds<0) ;
if ( isempty ( idx_s )==0)
alpha = min ( alpha ,min(-s( idx_s ) ./ ds ( idx_s ) ) ) ;
end
%Update x , z , s
x = x + eta * alpha*dx ;
z = z + eta * alpha*dz ;
s = s + eta * alpha*ds ;
k = k+1;
%Update rhs and mu
rL = G*x + g - A*z ;
rs = s - A' * x + b ;
rsz = s .* z ;
mu = sum( z .* s ) /nA;
end
xstop = x ;
zstop = z ;
sstop = s ;