function [ x,mu] = range_space(L,A, g ,b)
 % RANGE SPACE sol v e s the e qual i ty constrained convex QP:
 % min 1/2x 'Gx+g ' x (G i s required to be postiv e
% d e f i n i t e )
 % s . t . A' x = b (A i s required to have f u l l column
% rank )
 % where the number of v a r i ab l e s i s n and the number of co nst r ain ts i s m.
 % The range space of the OP i s used to fi nd the so l ut i on .
 %
 % Call
 % [ x , u ] = range space (L ,A, g , b)
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
 % mu : the lagrangian m u l t i p l i e r s

 % By : Carsten V\ ¨olcker , s961572 & Esben Lundsager Hansen , s022022 .
 % Subject : Numerical Methods for Sequential Quadratic Optimization ,
 % Master Thesis , IMM, DTU, DK-2800 Lyngby .
 % Supervisor : John Bagterp Jørgensen , Assistant P rofessor & Per Grove
% Thomsen , Professor .
 % Date : 08. february 2007.
 % Reference : ---------------------

 Lt = L' ;
 K = Lt\A;
 H = K'*K;
 w = Lt\g ;
 z = b+K'*w;
 M = chol (H) ;
 mu = M' \ z ;
 mu = M\mu;
 y = K*mu-w;
 x = L\y ;
end