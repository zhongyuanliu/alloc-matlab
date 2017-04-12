function [ x ,mu, info , perf ] = primal_active_set_method(G, g ,A, b , x , w_non , pbc , opts ,trace )

 % PRIMAL ACTIVE SET METHOD Solving an i n eq u al it y constrained QP of the
 % form :
 % min f(x) = 0.5* x'*G*x + g*x
 % s . t . A*x >= b ,
 % by sol vi ng a sequence of eq ual it y constrained QP's using the primal
 % a cti v e set method . The method uses the range space procedure or the nul l
 % space procedure to solv e the KKT system . Both the range space and the
 % n ull space procedures has been provided with f a c t o r i z a t i o n updates .
 %
 % Call
 % x = primal acti ve set method(G, g , A, b , w_non , pbc )
 % x = primal acti ve set method(G, g , A, b , w_non , pbc , opts )
 % [ x , mu, info , perf ] = pri mal acti ve set method( . . . )
 %
 % Input parameters
 % G : The Hessian matrix of the objective function , size nxn .
 % g : The l i n e a r term of the o bj ec t i v e function , size nx1 .
 % A : The co nstr ain t matrix holding the constraints , size nxm.
 % b : The right -hand side of the const raints , si z e mx1.
 % x : Starting point , size nx1 .
 % w_non : L ist of i n ac t i v e const rai nts , pointing on co nst r ai nt s in A.
 % pbc : L ist of corresponding constraints , pointing on c on str a i nt s in
 % A. Can be empty .
 % opts : Vector with 3 elements :
 % opts(1) = Tolerance used to s t a b i l i z e the methods numerically .
 % I f | value | <= opts(1) , then value i s regarded as zero .
 % opts(2) = maximum no . of it e r a t i o n steps .
 % opts(3) = 1 : Using null space procedure .
 % 2 : Using null space procedure with f a c t o r i z a t i o n
 % update .
 % 3 : Using null space procedure with f a c t o r i z a t i o n
 % update based on fi xe d and f r e e v ar i a bl e s . Can only
 % be called , if the i ne q ua l i ty constrained QP i s
 % setup on the form seen in QP solver .
 % I f opts i s not given or empty , the defa ult opts = [1 e-8 1000 3 ] .
 %
 % Output parameters
 % x : The optimal sol ut i on .
 % mu : The Lagrange m u l t i p l i e r s at the optimal sol ut i on .
 % i nfo : Performace information , vector with 3 elements :
 % i nf o(1) = f i n a l values of the ob je ct i v e function .
 % i nf o(2) = no . of it e r a t i o n steps .
 % i nf o(3) = 1 : Feasible so l ut i on found .
 % 2 : No. of it e r a t i o n steps exceeded .
 % perf : Performace , st r uct holding :
 % perf . x : Values of x , size i s nx( it +1) .
 % perf . f : Values of the o bj ec t i v e function , size i s 1x( it +1) .
 % perf .mu : Values of mu, size i s nx( it +1) .
 % perf . c : Values of c(x ) , size i s mx( it +1) .
 % perf .Wa : Active set , size i s mx( it +1) .
 % perf .Wi : Inac ti ve set , size i s mx( it +1) .
 %
 % By : Carsten V\ ¨olcker , s961572 .
 % Esben Lundsager Hansen , s022022 .
 % Subject : Numerical Methods f or Sequential Quadratic Optimization .
 % M. Sc . , IMM, DTU, DK-2800 Lyngby .
 % Supervisor : John Bagterp Jørgensen , Assistant Pr ofe ssor .
 % Per Grove Thomsen , Professor .
 % Date : 07. June 2007.

 % the size of the con st ra int matrix , where the c on st ra i nt s are given columnwise

 [ n ,m] = size(A) ;

 nb = 2*n ; % number of bounds
 ngc = m - nb ; % number of general co nst r ai nt s

 % i n it i a l i z e . . .
 z = zeros(m, 1 ) ;
 x0 = x ;
 f0 = objective(G, g , x0 ) ;
 mu0 = z ;
 c0 = constraints(A, b , x0 ) ;
 w_act0 = z ;
 w_non0 =( 1 : 1 :m)';
% i n it i a l i z e options . . .
 tol = opts(1) ;
 it_max = opts(2) ;
 method = opts(3) ;
 % i n it i a l i z e contai ners . . .
%  trace =( nargout > 3) ;
 perf = [ ] ;
 if trace
 X = repmat( zeros(n , 1 ) ,1 , it_max ) ;
 F = zeros(1 , it_max ) ;
 Mu = repmat( z , 1 , it_max ) ;
 C = repmat( z , 1 , it_max ) ;
 W_act = repmat( z , 1 , it_max ) ;
 W_non = repmat( z , 1 , it_max ) ;
 end

 % i n it i a l i z e counters . . .
 it = 0;
 Q = [ ] ; T = [ ] ; L = [ ] ; rem = [ ] ; % both f or null space update and fo r
% null space update FXFR

if method == 3 % nul l space with FXFR-update
 nb1 = n *2; % number of bounds
 else
 nb1 = 0;
 end

 nab = 0; % number of act i ve bouns
 P = eye(n) ;

 w_act = 3 ;%!!!!!!!!!!!!!!!!!added
 Q_old = [ ] ; T_old = [ ] ; L_old = [ ] ; 
 rem = [ ] ;

 % it e r a t e . . .
 stop = 0;
 while ~stop
 it = it + 1;
 if it >= it_max
 stop = 2; % maximum no it e r a t i o n s exceeded
 end

 % c a l l range/ nul l space procedure . . .
 mu = z ;

 if method == 1
 [ p , mu_] = null_space(G,A( : , w_act ) ,G*x+g , zeros( length( w_act ) ,1) ) ;
 end

 if method == 2
 [ p , mu_,Q,T, L] = null_space_update(G,A( : , w_act ) ,G*x+g , zeros( length( w_act) ,1) ,Q,T, L , rem ) ;
 end

 if method == 3
 Cr = G*x+g ;
 A = A( : , w_act( nab+1:end ) ) ;
 % ajust C( : , r ) to make it correspond to the f a c t o r i z a t i o n s of
% the
% Fixed v a ri ab l es( whenever -1 appears at v ar i abl e i C( : , r ) i
% should change sign )
 if nab % some bounds are in the act iv e se t
 u_idx = find( w_act > nb1/2 & w_act< nb1+1) ;
 var = n-nab+u_idx ;
 Cr( var ) = -Cr( var ) ;
 A_( var , : ) = -A_( var , : ) ;
 end
 [ p , mu_,Q,T, L] = null_space_updateFRFX(Q,T, L,G, A , Cr , zeros( length( w_act ),1) ,nab , rem-nab ) ;
 end

 mu( w_act ) = mu_;

 if norm(p) < tol
 if mu > -tol
 stop = 1; % sol ut i on found
 else
 % compute index j of bound/ co nst ra int to be removed . . .
 [dummy, rem ] = min(mu_) ;
 [ w_act, w_non, A, P, x, nab, G, g ] = remove_constraint(rem ,A, w_act , w_non , x, P, nb1 , nab , n ,G, g , pbc , b) ;
 end
 else
 % compute step length and index j of bound/ const rai nt to be appended . . .
 [ alpha , j ] = step_length(A, b , x , p , w_non , nb , n , tol ) ;
 if alpha < 1
 % make constrained step . . .
 x = x + alpha*p ;
 [ w_act, w_non, A, P, x, nab, G, g, Q] = append_constraint( j ,A, w_act , w_non , x,P, nb1 , nab , n ,G, g ,Q, pbc , b) ; % r i s index of A
 else
 % make f u l l step . . .
 x = x + p ;
 end
 end
 % c o l l e c t i n g output in containers . . .
 if trace
 if nb1 % method 3 i s used
 X( : , it ) = P'* x ;
 else
 X( : , it ) = x ;
 end
 F( it ) = objective(G, g , x ) ;
 Mu( : , it ) = mu;
 C( : , it ) = constraints(A, b , x ) ;
 W_act( w_act , it ) = w_act ;
 W_non( w_non , it ) = w_non ;
 end
 end

 if nb1 % method 3 i s used
 x = P'* x ;
 end

 % building i nf o . . .
 info = [ objective(G, g , x ) it stop ] ;

 % bui lding perf . . .
 if trace
 X = X( : , 1 : it ) ; X = [ x0 X] ;
 F = F( 1 : it ) ; F = [ f0 F ] ;
 Mu = Mu( : , 1 : it ) ; Mu = [mu0 Mu] ;
 C = C( : , 1 : it ) ; C = [ c0 C] ;
 W_act = W_act( : , 1 : it ) ; W_act = [ w_act0 W_act ] ;
 W_non = W_non( : , 1 : it ) ; W_non = [ w_non0 W_non ] ;
 perf = struct('x',{X} ,'f',{F} ,'mu',{Mu} ,'c',{C} ,'Wa',{ W_act } ,'Wi',{W_non}) ;
 end

 function [ alpha , j ] = step_length(A, b , x , p , w_non , nb , n , tol )
 alpha = 1; j = [ ] ;
 for app = w_non
 if app > nb
 fv = 1 : 1 : n ; % general const rai nt
 else
 fv = mod(app-1,n) +1; % index of fi x ed v ar i abe l
 end
 ap = A( fv , app )'*p( fv ) ;
 if ap < -tol
 temp =(b( app ) - A( fv , app)'*x( fv ) ) /ap ;
 if -tol < temp && temp < alpha
 alpha = temp ; % sm all est step length
 j = app ; % index j of bound to be appended
 end
 end
 end

 % function [ w_act , w_non ] = append constraint(b , w_act , w_non , j , pbc)
 % w_act = [ w_act j ] ; % append const rai nt j to ac ti ve set
 % w_non = w_non( find( w_non ~= j ) ) ; % remove c onstr ai nt j from nonactive set
 % if ~isinf(b( pbc( j ) ) )
 % w_non = w_non( find( w_non ~= pbc( j ) ) ) ; % remove c onstr ai nt pbc( j ) from
% nonactive set , if not unbounded
 % end

 % function [ w_act , w_non ] = remove constraint(b , w_act , w_non , j , pbc)
 % w_act = w_act( find( w_act ~= j ) ) ; % remove c onstr ai nt j from a cti v e set
 % w_non = [ w_non j ] ; % append const rai nt j to nonactive se t
 % if ~isinf(b( pbc( j ) ) )
 % w_non = [ w_non pbc( j ) ] ; % append con st rai nt pbc( j ) to nonactive set , if
% not unbounded
 % end

 function [ w_act, w_non, C, P, x, nab, G, g ] = remove_constraint( wi ,C, w_act , w_non , x ,P,nb , nab , n ,G, g , pbc , b) % wi i s index of w_act
 j = w_act( wi ) ;

 if j < nb+1 % j i s a bound and we have to
% r eorganiz e the v ar i a bl e s
 var1 = n-nab+1;
 var2 = n-nab+wi ;

 temp = C( var1 , : ) ;
 C( var1 , : ) = C( var2 , : ) ;
 C( var2 , : ) = temp ;

 temp = x( var1 ) ;
 x( var1 ) = x( var2 ) ;
 x( var2 ) = temp ;

 temp = P( var1 , : ) ;
 P( var1 , : ) = P( var2 , : ) ;
 P( var2 , : ) = temp ;

 temp = G( var1 , var1 ) ;
 G( var1 , var1 ) = G( var2 , var2 ) ;
 G( var2 , var2 ) = temp ;

 temp = g( var1 ) ;
 g( var1 ) = g( var2 ) ;
 g( var2 ) = temp ;
 nab = nab - 1;

 temp = w_act( wi ) ;
 w_act( wi ) = w_act(1) ;
 w_act(1) = temp ;
 j = w_act(1) ;
 end
 w_act = w_act( find( w_act ~= j ) ) ; % bound/ general const rai nt j i s
removed from a cti v e set
 w_non = [ w_non j ] ; % bound/ general const rai nt j
appended to nonactive set

 if ~isinf(b( pbc( j ) ) )
 w_non = [ w_non pbc( j ) ] ; % append bound/ con st rai nt pbc( j )
% to nonactive set , if not unbounded
 end

 function [ w_act, w_non, C, P, x, nab, G, g, Q] = append_constraint( j ,C, w_act , w_non , x ,P,nb , nab , n ,G, g ,Q, pbc , b) % j i s index of C
 if j < nb+1 % j i s a bound and we have to
% reorgani ze the v ar i a bl es
 var1 = find( abs(C( : , j ) )==1);
 var2 = n-nab ;

 temp = C( var1 , : ) ;
 C( var1 , : ) = C( var2 , : ) ;
 C( var2 , : ) = temp ;

 temp = Q( var1 , : ) ;
 Q( var1 , : ) = Q( var2 , : ) ;
 Q( var2 , : ) = temp ;

 temp = x( var1 ) ;
 x( var1 ) = x( var2 ) ;
 x( var2 ) = temp ;

 temp = P( var1 , : ) ;
 P( var1 , : ) = P( var2 , : ) ;
 P( var2 , : ) = temp ;

 temp = G( var1 , var1 ) ;
 G( var1 , var1 ) = G( var2 , var2 ) ;
 G( var2 , var2 ) = temp ;

 temp = g( var1 ) ;
 g( var1 ) = g( var2 ) ;
 g( var2 ) = temp ;
 nab = nab + 1;
 w_act = [ j w_act ] ; % j( i s a bound ) i s appended to
% ac ti ve set
 else
 w_act = [ w_act j ] ; % j( i s a general c onstr ai nt )
% i s appended to act iv e set
 end
 w_non = w_non( find( w_non ~= j ) ) ; % bound/ general c onstr ai nt j
% i s removed fom nonactive set
 if ~isinf(b( pbc( j ) ) )
 w_non = w_non( find( w_non ~= pbc( j ) ) ) ; % remove bound/ const rai nt pbc( j
% ) from nonactive set , if not unbounded
 end

 function f = objective(G, g , x )
 f = 0.5* x'*G*x + g'* x ;
 

 function c = constraints(A, b , x )
 c = A'* x - b ;

 function l = lagrangian(G, g ,A, b , x ,mu)
 L = objective(G, g ,A, b , x ,mu) - mu( : )'* constraints(G, g ,A, b , x ,mu) ;