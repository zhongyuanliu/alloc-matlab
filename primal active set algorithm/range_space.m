function [ x,lambda] = range_space(G,g,A,b)
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