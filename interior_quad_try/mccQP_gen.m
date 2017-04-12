function [ x_stop , y_stop , z_stop , s_stop , k ] = mccQP_gen( x , y , z , s ,G, g ,C, d ,A, b)
%??????????????????????????????????????????????????????????
% mccQP gen .m
%
% This f u n c t i on s o l v e s a QP problem o f the form
% min 0 . 5 x ¡¯Gx + g ¡¯ x
% s . t . A¡¯ x = b
% C¡¯ x >= d
% us ing the Mul t ipl e Ce n t r a l i t y Co r r e c t i on s (MCC) method .
%
% Input :
% x : s t a r t i n g point f o r the ve c tor x (n x 1 ve c tor )
% y : s t a r t i n g point f o r the ve c tor x (nA x 1 ve c tor )
% z : s t a r t i n g point f o r the ve c tor x (nC x 1 ve c tor )
% s : s t a r t i n g point f o r the s lack?ve c tor x (nC x 1 ve c tor )
% G: He s s ian (n x n matrix )
% g : (n x 1 ve c tor )
% C: l e f t hand s i d e o f i n e q u a l i t y c o n s t r a i n t s (n x nC matrix )
% d : r i gh t hand s i d e o f i n e q u a l i t y c on s t r a i n t s (nC x 1 ve c tor )
% A: l e f t hand s i d e o f e qu a l i t y c o n s t r a i n t s (n x nA matrix )
% b : r i gh t hand s i d e o f e qu a l i t y c o n s t r a i n t s (nA x 1 ve c tor )
% where nC and nA ar e the numbers o f i n e q u a l i t y and e qu a l i t y
% c o n s t r a i n t s .
% Output :
% x_stop : s o l u t i o n x
% y_stop
% z s t op
% s s t op
% k : number o f i t e r a t i o n s used
%
% Thomas Reslow Kr¡§uth , s021898
%??????????????????????????????????????????????????????????
B min = 0 . 1 ; B max = 1 0 ; d alpha = 0 . 1 ; gamma = 0 . 1 ;
e ta = 0 . 9 9 9 9 5 ;
%r e s i d u a l s ar e computed
[mA,nA] = s i z e (A) ;
[mC, nC] = s i z e (C) ;
e = one s (nC, 1 ) ;
rL = Gx + g ? Ay ? Cz ;
rA = ?A¡¯  x + b ;
rC = ?C¡¯  x + s + d ;
r s z = s .  z ;
mu = sum( z .  s ) /nC;
%k : number o f i t e r a t i o n s , e p s i l on : t o l e r a n c e s
k = 0 ;
maxk = 200;
eps L = 1e?10; eps A = 1e?10; eps C = 1e?10; eps mu = 1e?10;
whi l e ( k<=maxk && norm( rL )>=eps L && norm( rA)>=eps A && norm( rC)>=eps C
&& abs (mu)>=eps mu )
%Solve system with a Newton?l i k e method/ Fa c t o r i z a r i on
l h s = [G,?A,?C;?A¡¯ , s p a l l o c (nA, nA, 0 ) , s p a l l o c (nA, nC, 0 ) ;?C¡¯ , s p a l l o c (nC,
nA, 0 ) , sp a r s e (?diag ( s . / z ) ) ] ;
[L,D,P] = l d l ( l h s ) ;
rhs = [?rL;?rA;?rC+r s z . / z ] ;
dxyz a = P(L¡¯ \ (D\(L\(P¡¯  rhs ) ) ) ) ;
dz a = dxyz a ( l ength ( x )+l ength ( y )+1: l ength ( x )+l ength ( y )+l ength ( z ) ) ;
ds a = ?(( r s z+s .  dz a ) . / z ) ;
%Compute a l p h a a f f
alpha a = 1 ;
i d x z = f i n d ( dz a <0) ;
i f ( isempty ( i d x z )==0)
alpha a = min ( alpha a ,min(?z ( i d x z ) . / dz a ( i d x z ) ) ) ;
end
i d x s = f i n d ( ds a <0) ;
i f ( isempty ( i d x s )==0)
alpha a = min ( alpha a ,min(?s ( i d x s ) . / ds a ( i d x s ) ) ) ;
end
%Compute the a f f i n e d u a l i t y gap
mu a = ( ( z+alpha a  dz a ) ¡¯( s+alpha a  ds a ) ) /nC;
%Compute the c e n t e r i n g parameter
sigma = (mu a/mu) ? 3 ;
mu t = sigma mu;
r s z = r s z + ds a .  dz a ? mu te ;
rhs = [?rL;?rA;?rC+r s z . / z ] ;
dxyz p = P(L¡¯ \ (D\(L\(P¡¯  rhs ) ) ) ) ;
dx p = dxyz p ( 1 : l ength ( x ) ) ;
dy p = dxyz p ( l ength ( x )+1: l ength ( x )+l ength ( y ) ) ;
dz p = dxyz p ( l ength ( x )+l ength ( y )+1: l ength ( x )+l ength ( y )+l ength ( z ) ) ;
ds p = ?(( r s z+s .  dz p ) . / z ) ;
%Compute alpha
alpha p = 1 ;
i d x z = f i n d ( dz p<0) ;
i f ( isempty ( i d x z )==0)
alpha p = min( alpha p ,min(?z ( i d x z ) . / dz p ( i d x z ) ) ) ;
end
i d x s = f i n d ( ds p <0) ;
i f ( isempty ( i d x s )==0)
alpha p = min( alpha p ,min(?s ( i d x s ) . / ds p ( i d x s ) ) ) ;
end
%Modi f ied c e n t e r i n g d i r e c t i o n s
alpha h = min( alpha p+d alpha , 1 ) ;
%Maximum number o f c o r r e c t i o n s
K=2;
%Co r r e c t i on s
j = 1 ;
whi l e ( j<=K)
%Compute t r i a l point
z t = z + alpha h dz p ;
s t = s + alpha h  ds p ;
%De f ine t a r g e t
v t h i l d e = s t .  z t ;
v tmp = v t h i l d e ;
f o r i =1: l ength ( v tmp )
i f ( v tmp ( i )<B minmu t )
v tmp ( i ) = B minmu t ;
e l s e i f ( v tmp ( i )>B maxmu t )
v tmp ( i ) = B maxmu t ;
end
end
v t = v tmp ;
%Compute c o r r e c t o r
r s z = ?( v t ? v t h i l d e ) ;
f o r i =1: l ength ( r s z )
i f ( r s z ( i )<(?B maxmu t ) )
r s z ( i ) = ?B maxmu t ;
end
end
rhs = [ z e r o s ( l ength ( rL ) ,1) ; z e r o s ( l ength ( rA) ,1) ; r s z . / z ] ;
dxyz m = P(L¡¯ \ (D\(L\(P¡¯  rhs ) ) ) ) ;
dx m = dxyz m ( 1 : l ength ( x ) ) ;
dy m = dxyz m( l ength ( x )+1: l ength ( x )+l ength ( y ) ) ;
dz m = dxyz m( l ength ( x )+l ength ( y )+1: l ength ( x )+l ength ( y )+l ength ( z
) ) ;
ds m = ?(( r s z+s .  dz m) . / z ) ;
%Composite d i r e c t i o n
dx = dx p + dx m;
ds = ds p + ds m;
dy = dy p + dy m;
dz = dz p + dz m;
%Compute alpha
alpha = 1 ;
i d x z = f i n d ( dz<0) ;
i f ( isempty ( i d x z )==0)
alpha = min ( alpha ,min(?z ( i d x z ) . / dz ( i d x z ) ) ) ;
end
i d x s = f i n d ( ds<0) ;
i f ( isempty ( i d x s )==0)
alpha = min ( alpha ,min(?s ( i d x s ) . / ds ( i d x s ) ) ) ;
end
%Test f o r improvement
i f ( alpha >= alpha p + gamma d alpha )
j=j +1;
dx p = dx ;
ds p = ds ;
dz p = dz ;
dy p = dy ;
alpha p = alpha ;
alpha h = min ( alpha p+d alpha , 1 ) ;
e l s e
dx = dx p ;
ds = ds p ;
dz = dz p ;
dy = dy p ;
j =10;
end
end %inne r whi le?loop
%Update x , y , z , s
x = x + e ta  alphadx ;
y = y + e ta  alphady ;
z = z + e ta  alphadz ;
s = s + e ta  alphads ;
k = k+1;
%Update rhs
rL = Gx + g ? Ay ? Cz ;
rA = ?A¡¯  x + b ;
rC = ?C¡¯  x + s + d ;
r s z = s .  z ;
mu = sum( z .  s ) /nC;
end
%Output
x_stop = x ;
y_stop = y ;
z s t op = z ;
s s t op = s ;