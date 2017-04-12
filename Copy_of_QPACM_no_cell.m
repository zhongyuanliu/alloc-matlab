function [x,flag]=QPACM(B,b,C,c,A,d,x0)
%%
% min f(x)=1/2*x'*B*x - b'*x  %NB
% C*x = c
% A * x <= d %NB
flag = 0;
tol = 1e-6;
x=x0;
% W.A = A([1],:);
% W.d = d([1]);
W.A=zeros(0,size(A,2));%empty
W.d=[];%
% [CW.A,IA] = setdiff(A,W.A,'rows','stable');%set1 - set2
CW.A = A;
CW.d = d;
% CW.d = d(IA);
m_C = size(C,1);
max_iter = 200;
i = 0;
while i<max_iter
i = i + 1;
flag = i;
% check_rank = (rank([C;W.A])==size([C;W.A],1));
[ p,r] = null_space(B,B*x-b,[C;W.A],0*[c;W.d]);%null space is better
p = p .* (abs(p)>tol);
r = r .* (abs(r)>tol);
% lambda = r(1:m_C);
miu = r(m_C + 1:end);

if sum(abs(p)) == 0
    if min(miu)>=0
        flag=i;
        break;
    else
        [~,I] = min(miu);
        CW.A = [CW.A;W.A(I,:)];
        CW.d = [CW.d;W.d(I)];
        if (~isempty(W.A(I,:))) && (~isempty(W.d(I)))
        W.A(I,:) = [];%delete constr
        W.d(I) = [];
        end
    end
else
    temp = CW.A * p;
    I = find(temp<-tol);
    
    if isempty(I)
        alpha = 1;
    else
    [alpha,Imin_alpha] = min((CW.A(I,:) * x - CW.d(I))./(CW.A(I,:) * p));
    [alpha] = min(1,alpha);
    end
    x = x -alpha * p;%NB
    if  alpha<1    %alpha>0 &&
%         [~,Imin_alpha] = min((CW.A(I,:) * x - CW.d(I))./(CW.A(I,:) * p));
        W.A = [W.A;CW.A(I(Imin_alpha),:)];%add constr
        W.d = [W.d;CW.d(I(Imin_alpha))];
        CW.A(I(Imin_alpha),:) = [];%delete
        CW.d(I(Imin_alpha)) = [];
    end
end
end
% 1/2*x'*B*x - b'*x;
end
function [ x,lambda] = range_space(G,g,A,b)
%%
 % RANGE SPACE sol v e s the e qual i ty constrained convex QP:
 % min 1/2*x'Gx+g ' x (G i s required to be postiv e
% d e f i n i t e )
 % s . t . A x = b (A i s required to have f u l l column
% rank )
% G*x +g +A*lambda
 % where the number of v a r i ab l e s i s n and the number of co nst r ain ts i s m.
 % The range space of the OP i s used to fi nd the so l ut i on .
 % Input parameters
 % L : i s the cholesky f a c t o r i z a t i o n of the Hessian matrix (
% nxn ) of the QP.
 % A : i s the const rai nt matrix (nxm) : every column contains
% a from the
 % e qua li ty : a ' x = b .
 % g : i s the gradient ( nx1 ) of the QP.
 % b : i s the ri ght hand si de of the co nst r ai nt s .
 % Output parameters
 % x : the sol u ti o n
 % lambda : the lagrangian m u l t i p l i e r s
 % Date : 19. July 2016.
 % Reference : Numerical Algorithms for Sequential Quadratic
 % Optimization,pange 10
%  A = -A;
 g=-g;
%  b = -b;
if isempty(A)
     x = -G\g ;
     lambda = [];
else
 L =chol(G,'lower');%L*L' = G
 K = (L)\A';
 H = K'*K;
 w = L\g ;
 z = b+K'*w;
 M = chol (H,'lower');
 q = M \ z ;
 lambda = M'\q;
 r = K*lambda-w;
 x = L'\r ;
%  lambda = -lambda;
end
end

function [ x , u ]= null_space(G,g,A,b)
%%
% NULL SPACE sol v e s the eq ual i ty constrained convex QP:
% min 1/2x 'Gx+g ' x (G i s required to be postiv e semi
% d e f i n i t e )
% s . t . A' x = b (A i s required to have f u l l column
% rank )
% where the number of v a ri a b l es i s n and the number of c o nst ra int s i s m.
% The nul l space of the OP i s used to fi nd the so l ut io n .
%
% Call
% [ x , u ] = n ul l sp ac e (G,A, g , b)
%
% Input parameters
% G : i s the Hessian matrix (nxn ) of the QP.
% A : i s the c onstr ai nt matrix (nxm) : every column contains
% a from the
% e qual i ty : a ' x = b .
% g : i s the gradient (nx1 ) of the QP.
% b : i s the r i ght hand si de of the c o nst ra i nt s .
%
% Output parameters
% x : the sol u ti o n
% mu : the lagrangian m u l t i p l i e r s
%
% By : Carsten V\ ¨olcker , s961572 & Esben Lundsager Hansen , s022022 .
% Subject : Numerical Methods f or Sequential Quadratic Optimization ,
% Master Thesis , IMM, DTU, DK-2800 Lyngby .
% Supervisor : John Bagterp Jørgensen , Assistant Pr ofe ssor & Per Grove
% Thomsen , Professor .
% Date : 08. february 2007.
g=-g;
A=A';
[ n ,m] = size (A) ;
if ( m~=0 ) % for si t u a t i o ns where A i s empty
    [Q,R] = qr (A) ;
    Q1 = Q( : , 1 :m) ;
    Q2 = Q( : ,m+1:n) ;
    R = R( 1:m, : ) ;
    py = R'\b ;
    Q2t = Q2' ;
    gz = Q2t*(G*(Q1*py ) + g ) ;
    Gz = Q2t*G*Q2;
    L = chol(Gz)' ;
    pz = L\-gz ;
    pz = L'\ pz ;
    x = Q1*py + Q2*pz ;
    u = R\(Q1' * (G*x + g ) ) ;
    else
    x = -G\g ;
    u = [ ] ;
end
end