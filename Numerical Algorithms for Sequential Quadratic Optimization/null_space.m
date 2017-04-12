function [ x , u ]= null_space(G,A,g,b)

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