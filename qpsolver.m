function [solution,alloc_out,thruster_data,status]=qpsolver(init,thruster_data,T_r,N_enabled_thruster,rudder_table0,no_azi_angle_constr,method)
% method == 1 qpkwik
% method == 2 sef made primal active set

% solves a quadratic program to determine an optimal solution, x.
% It minimizes the quadratic objective function, J = 0.5*x'*H*x + f'*x,
% subject to linear inequality constraints, A*x >= b, and linear equality
% constraints, Aeq*x = beq, where x is a column vector of length n.
%
% H is the n-by-n hessian matrix, which must be symmetric and positive
% definite.  L is its lower-triangular Cholesky decomposition and Linv is
% the inverse of L.  Consequently, Linv can be computed from H as follows
% in MATLAB:
%   L = chol(H, 'lower')
%   Linv = L\eye(size(H,1))
% Note that H = L*L'
%
% Input arguments (all are mandatory):
%    Linv:  n-by-n (n>0) matrix representing the inverse of L.
%       f:  column vector of length n (n>0).
%       A:  m-by-n matrix of linear inequality constraint coefficients.  If
%           your problem has no inequality constraints, use [].
%       b:  column vector of length m, the right-hand side of A*x >= b.  If
%           your problem has no inequality constraints, use "zeros(0,1)".
%     Aeq:  q-by-n matrix of linear equality constraint coefficients,
%           q <= n.  If your problem has no equality constraints, use [].
%           NOTE:  equality constraints must be linearly independent with
%           rank(Aeq) = q.
%     beq:  column vector of length q, the right-hand side of Aeq*x = beq.
%           If your problem has no equality constraints, use "zeros(0,1)".
%     iA0:  logical vector of length m.  If your problem has no inequality
%           constraints, use "false(0,1)".  For a "cold start": use
%           false(m,1).  For a "warm start": if iA0(i) == true, the
%           algorithm begins with A(i,:)*x = b(i), i.e., with the ith
%           inequality active.
%           Tips:
%           (1) Normally, you should use the optional output argument iA
%           from a previous solution as the input iA0 in the next
%           calculation.
%           (2) if iA0(i) and iA0(j) are true, rows i and j of A should be
%           linearly independent.  Otherwise, the solution may fail (status
%           = -2). This should not happen if you use iA from a previous
%           solution as recommended above.
% options:  Structure with following fields used by the QP solver.
%
%                  DataType: string, either 'double' or 'single'.  All the
%                            real input arguments to the "mpcqpsolver"
%                            command must match this data type and it is
%                            used in both simulation and code generation.
%                            Default value is 'double'.
%                   MaxIter: scalar, maximum number of iterations allowed
%                            when computing QP solution.  Default value is
%                            200.
%            FeasibilityTol: scalar, tolerance used to verify that
%                            inequality constraints are satisfied at the
%                            optimal solution.  Increasing this causes
%                            MPCQPSOLVER to allow larger constraint
%                            violations. Default value is 1.0e-6.
%
%           Use the MPCQPSOLVEROPTIONS command to create such a structure
%           with default values.
%
% Output arguments:
%       x: column vector of length n, representing the optimal solution.
%          MPCQPSOLVER will return x in all cases, but it is likely to be
%          sub-optimal and/or infeasible unless status > 0.
%  status: scalar, indicating the validity of the returned x as follows:
%            >  0:  x is optimal, representing the number of iterations
%                   used in optimization.
%           ==  0:  x was obtained when maximum number of iterations is
%                   reached.  It may be sub-optimal and/or violate A*x >=b.
%           == -1:  The problem appears to be infeasible.  In other words,
%                   A*x >= b cannot be satisfied.
%           == -2:  An unrecoverable numerical error occurred (e.g., see
%                   description of iA0).
%      iA: logical vector of length m, indicating the inequalities that
%          are active (at equality) at the solution.  In a series of
%          problems in which A is constant, you can use iA from one
%          solution as the input iA0 to the next ("warm start").
%  lambda: structure of Lagrange multipliers with two fields:
%           ineqlin: column vector of length m, multipliers of the
%                    inequality constraints.  Must be non-negative at an
%                    optimal solution.
%             eqlin: column vector of length q, multipliers of the equality
%                    constraints.  No sign restriction.

%  command for code generation:
%  cfg = coder.config('lib');
%  codegen -config cfg -c  -args {init,thruster_data,T_r,N_enabled_thruster,rudder_table,0,1} qpsolver

%  update:
%  1.01 add check_phi in pre_qp.m to make sure initial phi meets
%  constraint. |29-11-2016 by Liu
%  1.02 correct mark8 and add if status < 0 then break to avoid risk of
%  infinite loop. |29-11-2016 by Liu
%  1.04 change alloc_out(i).phi = thruster_data(i).phi;%output is equal to input for tunnel thruster
%  1.05 improve function check_angle(phi_min,phi_,phi,phiPlus,phi_max) and angleMaxMin(angle_start,angle_end,phi,phi1,phi2)
%  1.06 guarantee initial x feasible, add mark8(i) = 2 to avoid infinite
%  loop, add count_run > 50 then break |03-01-2017 by Liu
%  1.07 guarantee tunnel angle input feasible |23-01-2017 by Liu
%  4.01 base thrust added, Tmin of azi must be zero in thrust_config !!!!!
%  and bastT must be half of real Tmin. dt*dphi_max is replaced by a bigger constant |01-02-2017 by Liu
%  4.04 azi use efficiency one by one |03-02-2017 by Liu
%  4.05 azi efficiency is smoothy |06-02-2017 by Liu
%  4.06 u is scaled. c_scale * Tmax^2 == 1 |13-03-2017 by Liu
coder.inline('never');
rudder_table = rudder_table0';%!!!!!!!!!!!!!!!
% no_azi_angle_constr = 1;
%% rudder table processing
persistent Rangle_T Rangle_Fangle Fangle_Rangle rudder_dat_init;
persistent on_border phi_guider_old;
% persistent iA;
if isempty(rudder_dat_init)%assume this is only one type of rudder
    rudder_dat_init = false;
    N_poly = 2;
    [Rangle_T, ~] = polyfit(rudder_table(:,1), rudder_table(:,5),N_poly);%rudder angle vs T
    [Rangle_Fangle, ~] = polyfit(rudder_table(:,1), rudder_table(:,4),N_poly);%rudder angle vs F angle
    [Fangle_Rangle, ~] = polyfit(rudder_table(:,4), rudder_table(:,1),N_poly);%F angle vs rudder angle
end

if isempty(on_border)
    on_border = zeros(8,1);
    phi_guider_old = zeros(8,1);
end
max_semi_angle = 15/180*pi;%max half range of phi;
NA_fpp = size(rudder_table,1) * 2;
% coder.varsize('alloc_out');%!!!!!!!!!very important???
alloc_out = struct ('label',{0,0,0,0,0,0,0,0}, ...
    'enable',0, ...
    'TSP',0, ...
    'ASP',0, ...
    'Tx',0, ...
    'Ty',0, ...
    'T',0, ...
    'phi',0, ...
    'Tm',0);

Ce0=[0,0.800000000000000;0.174532925199433,0.820000000000000;0.349065850398866,0.840000000000000;0.523598775598299,0.870000000000000;0.698131700797732,0.900000000000000;0.872664625997165,0.910000000000000;1.04719755119660,0.930000000000000;1.22173047639603,0.960000000000000;1.39626340159546,0.750000000000000;1.57079632679490,0.500000000000000;1.74532925199433,0.750000000000000;1.91986217719376,1;2.09439510239320,1;2.26892802759263,1;2.44346095279206,1;2.61799387799149,1;2.79252680319093,1;2.96705972839036,1;3.14159265358979,1;3.31612557878923,1;3.49065850398866,1;3.66519142918809,1;3.83972435438753,1;4.01425727958696,1;4.18879020478639,1;4.36332312998582,1;4.53785605518526,1;4.71238898038469,1;4.88692190558412,0.980000000000000;5.06145483078356,0.960000000000000;5.23598775598299,0.930000000000000;5.41052068118242,0.910000000000000;5.58505360638185,0.890000000000000;5.75958653158129,0.870000000000000;5.93411945678072,0.840000000000000;6.10865238198015,0.820000000000000;0,0.800000000000000];
Ce3=ones(size(Ce0));
Ce3(:,1)=Ce0(:,1);
Ce2 = [0,1;0.174532925199433,1;0.349065850398866,1;0.523598775598299,1;0.698131700797732,1;0.872664625997165,1;1.04719755119660,1;1.22173047639603,1;1.39626340159546,1;1.57079632679490,1;1.74532925199433,0.980000000000000;1.91986217719376,0.960000000000000;2.09439510239320,0.930000000000000;2.26892802759263,0.910000000000000;2.44346095279206,0.890000000000000;2.61799387799149,0.870000000000000;2.79252680319093,0.840000000000000;2.96705972839036,0.820000000000000;3.14159265358979,0.800000000000000;3.31612557878923,0.820000000000000;3.49065850398866,0.840000000000000;3.66519142918809,0.870000000000000;3.83972435438753,0.900000000000000;4.01425727958696,0.910000000000000;4.18879020478639,0.930000000000000;4.36332312998582,0.960000000000000;4.53785605518526,0.750000000000000;4.71238898038469,0.500000000000000;4.88692190558412,0.750000000000000;5.06145483078356,1;5.23598775598299,1;5.41052068118242,1;5.58505360638185,1;5.75958653158129,1;5.93411945678072,1;6.10865238198015,1;0,1];
Ce1 = [0,1;0.174532925199436,1;0.349065850398866,1;0.523598775598297,1;0.698131700797736,1;0.872664625997166,1;1.04719755119660,1;1.22173047639603,1;1.39626340159547,0.750000000000000;1.57079632679490,0.500000000000000;1.74532925199433,0.750000000000000;1.91986217719377,0.960000000000000;2.09439510239320,0.930000000000000;2.26892802759263,0.910000000000000;2.44346095279206,0.900000000000000;2.61799387799150,0.870000000000000;2.79252680319093,0.840000000000000;2.96705972839036,0.820000000000000;3.14159265358980,0.800000000000000;3.31612557878923,0.820000000000000;3.49065850398866,0.840000000000000;3.66519142918810,0.870000000000000;3.83972435438753,0.890000000000000;4.01425727958696,0.910000000000000;4.18879020478639,0.930000000000000;4.36332312998583,0.960000000000000;4.53785605518526,0.980000000000000;4.71238898038469,1;4.88692190558413,1;5.06145483078356,1;5.23598775598299,1;5.41052068118242,1;5.58505360638185,1;5.75958653158129,1;5.93411945678072,1;6.10865238198015,1;0,1];
Ce1 = [0,0.800000000000000;0.174532925199433,0.820000000000000;0.349065850398866,0.840000000000000;0.523598775598299,0.870000000000000;0.698131700797732,0.900000000000000;0.872664625997165,0.910000000000000;1.04719755119660,0.930000000000000;1.22173047639603,0.960000000000000;1.39626340159546,0.750000000000000;1.57079632679490,0.500000000000000;1.74532925199433,0.750000000000000;1.91986217719376,1;2.09439510239320,1;2.26892802759263,1;2.44346095279206,1;2.61799387799149,1;2.79252680319093,1;2.96705972839036,1;3.14159265358979,1;3.31612557878923,1;3.49065850398866,1;3.66519142918809,1;3.83972435438753,1;4.01425727958696,1;4.18879020478639,1;4.36332312998582,1;4.53785605518526,1;4.71238898038469,1;4.88692190558412,0.980000000000000;5.06145483078356,0.960000000000000;5.23598775598299,0.930000000000000;5.41052068118242,0.910000000000000;5.58505360638185,0.890000000000000;5.75958653158129,0.870000000000000;5.93411945678072,0.840000000000000;6.10865238198015,0.820000000000000;0,0.800000000000000];
% Ce2 = Ce1;
% Ce3 = Ce1;
Ce1 = Ce3;
Ce2 = Ce3;
Ce = struct('label',{1,2,3,4,5,6,7,8},...
    'reserve',0, ...
    'value',{Ce1,Ce2,Ce3,Ce3,Ce3,Ce3,Ce3,Ce3});


dt = thruster_data(1).dt;
solution = zeros(N_enabled_thruster * 2 + 3,1);
phi_0 = zeros(8,1);
% phi_reserve = zeros(8,1);
% T_reserve = zeros(8,1);
phiPlus0 = zeros(8,1);
inside_fz = zeros(8,1);
alloc_out_T = zeros(8,1);
alloc_out_phi = zeros(8,1);

tol = 1e-6;



diff_phi = 1;
for step = 1:2 %step 1: global solution, step 2: local solution
    
    for i = 1:N_enabled_thruster %added v1.09
        % 0 no angle constraint no thrust constraint
        % 2 only angle constraint
        % 1 only thrust constraint
        % 3 means angle must be fixed because solution is inside fb zone
        if step == 1
            thruster_data(i).constr = 0;%no constraint at all
        else
            thruster_data(i).constr = 2;%no constraint at all
        end
    end
    
    status = 0;
    counter = 0;
    phi0 = zeros(8,1);
    run = true;

    while run
        counter = counter + 1;
        flag_run = 1;
        
        count_run = 0;
        
        while flag_run
            flag_run = 0;

            count_run = count_run + 1;
            [thruster_data,Nvar,H,f,A,b,Aeq,beq,x0] = pre_qp(init,thruster_data,T_r,N_enabled_thruster,Ce,NA_fpp,rudder_table,Rangle_T,Rangle_Fangle,no_azi_angle_constr,method);
            init = false;
            MaxIter=100;
            FeasibilityTol=1e-8;
            
            [L,p] = Chol_fc(H);
            % before chol
            if p==0
                switch method
                    case 1
                        Linv = inv(L);
                        opt = mpcqpsolverOptions;
                        opt.MaxIter = MaxIter;
                        opt.IntegrityChecks = false;
                        opt.FeasibilityTol = FeasibilityTol;
                        iA0 = false(size(b));
                        %                         if isempty(iA)
                        %                             iA = iA0;
                        %                         end
                        [solution,status,~] = Mpcqpsolver(Linv,f,A,b,Aeq,beq,iA0,opt);%faster to use iA instead of iA0
                        %     [solution1,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq)
                    case 2
                        %     x0=[x0;0;0;0];
                        col_Aeq=size(Aeq,2);
                        %     ????????????????ux,uy,?????????s1?s2?s3???
                        x0(col_Aeq - 2 :end) = Aeq(1:3,col_Aeq - 2 :end)\(beq(1:3) - Aeq(1:3,1:col_Aeq - 3) * x0(1:col_Aeq - 3));%does not work if there is tunnel!!!!!!!!!
                        
                        [solution,status] = QPACM(H,-f,Aeq,beq,-A,-b,x0);
                        % x = quadprog(H,f,A,b,Aeq,beq);
                end
            else
                solution = zeros(Nvar,1);
                status = -2;
            end
            
            %             solution = solution.*(abs(solution)>=FeasibilityTol);%!!!!!avoid nuemrical problem
            
            % output
            numOfTunnel = 0;
            numOfAzi = 0;
            numOfFpp = 0;
            count_target_thr = 0;
            for i=1:N_enabled_thruster
                alloc_out(i).label = thruster_data(i).label;
                alloc_out(i).enable = thruster_data(i).enable;
                
                switch thruster_data(i).type
                    case 4
                        temp = numOfTunnel * 2 + numOfAzi * 2 + numOfFpp * 2;
                        alloc_out(i).Tx = solution(temp + 1) * thruster_data(i).Tmax^2;
                        alloc_out(i).Ty = solution(temp + 2) * thruster_data(i).Tmax^2;
                        alloc_out(i).T = alloc_out(i).Ty;
                        %                 alloc_out(i).phi = atan2(solution(numOfTunnel * 2 + numOfAzi * 2 + numOfFpp * 2 + 2), ...
                        %                     solution(numOfTunnel * 2 + numOfAzi * 2 + numOfFpp * 2 + 1));
                        alloc_out(i).phi = thruster_data(i).phi;
                        alloc_out(i).phi = thruster_data(i).phi * (abs([alloc_out(i).T]) <= FeasibilityTol) ...
                            + alloc_out(i).phi .* (abs(alloc_out(i).T) > FeasibilityTol);
                        alloc_out(i).Tm = alloc_out(i).Ty * (thruster_data(i).x ...
                            - thruster_data(i). x0) + alloc_out(i).Tx ...
                            * (thruster_data(i).y0 - thruster_data(i).y);
                        
                        numOfTunnel = numOfTunnel + 1;
                    case 3
                        temp = numOfTunnel * 2 + numOfAzi * 2 + numOfFpp * 2;
                        alloc_out(i).Tx = solution(temp + 1) * thruster_data(i).Tmax^2;
                        alloc_out(i).Ty = solution(temp + 2) * thruster_data(i).Tmax^2;
                        alloc_out(i).T = sqrt(alloc_out(i).Tx^2 + alloc_out(i).Ty^2);
                        alloc_out(i).phi = atan2(alloc_out(i).Ty, ...
                            alloc_out(i).Tx);
                        alloc_out(i).phi = thruster_data(i).phi * (abs([alloc_out(i).T]) <= FeasibilityTol) ...
                            + alloc_out(i).phi .* (abs(alloc_out(i).T) > FeasibilityTol);
                        alloc_out(i).Tm = alloc_out(i).Ty * (thruster_data(i).x ...
                            - thruster_data(i). x0) + alloc_out(i).Tx ...
                            * (thruster_data(i).y0 - thruster_data(i).y);
                        
                       if step == 1 ...
                               && abs(thruster_data(i).phi_max - thruster_data(i).phi_min)> FeasibilityTol%means there is forbidden zone
                        if angle_dist(thruster_data(i).phi_min,alloc_out(i).phi) ...
                                > angle_dist(thruster_data(i).phi_min,thruster_data(i).phi_max) ...
                                + FeasibilityTol                  % if solution is inside fb zone
                            
%                             thruster_data(i).constr = 3;%means angle is fixed and  solution is inside fb zone
                            thruster_data(i).constr = 2;%added v4.00!!!!!!!!!!!!!!!
                            if on_border(i) == 0 %esle use previous phi, phi_ and phiPlus
                                if angle_dist(thruster_data(i).phi_max,alloc_out(i).phi) ...%phi closer to min
                                        > angle_dist(alloc_out(i).phi,thruster_data(i).phi_min)

                                    thruster_data(i).phi_ = thruster_data(i).phi_min;
                                    thruster_data(i).phiPlus = thruster_data(i).phi_min ...
                                        + 2 * dt * thruster_data(i).dphi_max;
                                    thruster_data(i).phi = thruster_data(i).phi_min ...
                                        + dt * thruster_data(i).dphi_max;
                                else %phi closer to max
                                    thruster_data(i).phi_ = thruster_data(i).phi_max ...
                                        - 2 * dt * thruster_data(i).dphi_max;
                                    thruster_data(i).phiPlus = thruster_data(i).phi_max;
                                    thruster_data(i).phi = thruster_data(i).phi_max ...
                                        - dt * thruster_data(i).dphi_max;
                                end                                
                            end                                                                
                            thruster_data(i).T_ = thruster_data(i).Tmin;
                            thruster_data(i).Tplus = thruster_data(i).Tmax;
                            
                            flag_run = 1;
                            
                        end
                        if on_border(i) == 1 
%                             if i == 1
%                                 temp = phi_guider_old(i);
%                                 save('phi_guider_old.txt','temp','-ASCII','-append');
%                             end
                            thruster_data(i).constr = 2;%added!!!!!
                            if abs(thruster_data(i).phi - phi_guider_old(i))>FeasibilityTol
                                thruster_data(i).phi = phi_guider_old(i);
                                thruster_data(i).phi_ = thruster_data(i).phi ...
                                    - 1e-4*0;
                                thruster_data(i).phiPlus = thruster_data(i).phi ...
                                    + 1e-4*0;

                                flag_run = 1;
                            else
                                thruster_data(i).phi_ = thruster_data(i).phi ...
                                - 1e-4*0;
                                thruster_data(i).phiPlus = thruster_data(i).phi ...
                                + 1e-4*0;
                            end
                            thruster_data(i).T_ = thruster_data(i).Tmin;
                            thruster_data(i).Tplus = thruster_data(i).Tmax;

                        end
                       end
                       
                        
                        
                        numOfAzi = numOfAzi + 1;
                        
                        
                    case 2
                        temp = numOfTunnel * 2 + numOfAzi * 2 + numOfFpp * 2;
                        alloc_out(i).Tx = solution(temp + 1) * thruster_data(i).Tmax^2;
                        alloc_out(i).Ty = solution(temp + 2) * thruster_data(i).Tmax^2;
                        alloc_out(i).T = sqrt(alloc_out(i).Tx^2 + alloc_out(i).Ty^2);
                        alloc_out(i).phi = atan2(alloc_out(i).Ty , ...
                            alloc_out(i).Tx);
                        alloc_out(i).phi = thruster_data(i).phi * (abs([alloc_out(i).T]) <= FeasibilityTol) ...
                            + alloc_out(i).phi .* (abs(alloc_out(i).T) > FeasibilityTol);
                        alloc_out(i).phi = sign(alloc_out(i).phi) * polyval(Fangle_Rangle,abs(alloc_out(i).phi));
                        alloc_out(i).Tm = alloc_out(i).Ty * (thruster_data(i).x ...
                            - thruster_data(i). x0) + alloc_out(i).Tx ...
                            * (thruster_data(i).y0 - thruster_data(i).y);
                        numOfFpp = numOfFpp + 1;
                end
            end
            if status < 0
                break;
            end
            
            if count_run > 50%avoid too many loops
                break;
            end
        end
        
        if step == 1
            
            diff_phi = 0;
            for i=1:N_enabled_thruster
                %         thruster_data(i).T = alloc_out(i).T;
                if thruster_data(i).type == 3
                    thruster_data(i).phi = alloc_out(i).phi;
                    diff_phi = diff_phi + mod(abs(phi0(i) - thruster_data(i).phi),pi);%check!!!!!!!!!
                    phi0(i) = thruster_data(i).phi;
                    
                end
            end
            
            error = (norm(solution(end-2:end)));
            if diff_phi < 0.01 || counter > 10 %this is the last loop in step one
                run = false;%break the while loop
                
                use_base = false;
                for i=1:N_enabled_thruster
                    if thruster_data(i).type == 3 ...
                            && alloc_out(i).T < thruster_data(i).base(1)/2%smaller than half base T
                        use_base = true;
                    end
                end
%                 use_base = true;%use base will increase error when base leads phi in fbzone
                if use_base
                    check_base_infb = 0;
                    for i=1:N_enabled_thruster
                        if thruster_data(i).type == 3
                            [alloc_out_T(i),alloc_out_phi(i),base_infb] = baseThrust(alloc_out(i).T,alloc_out(i).phi ...
                                ,thruster_data(i).base(1), thruster_data(i).base(2), ...
                                thruster_data(i).phi_min, thruster_data(i).phi_max);
                            check_base_infb = check_base_infb + base_infb;
                        end
                    end
                    % the drawback is that when in step
                    % 1 check_base_infb == 1, the result will
                    % not be accurate because phi
                    % in fbzone from step 1 will guide phi in step 2 to the
                    % fbzone border. If phi in step 2 is guided to fbzone,
                    % the result is accurate. However it is not allowed
                    if check_base_infb == 0 %only use base when it does not bring phi to fbzone
                        for i=1:N_enabled_thruster
                            if thruster_data(i).type == 3
                                alloc_out(i).T = alloc_out_T(i);
                                alloc_out(i).phi = alloc_out_phi(i);
                            end
                        end
                    end
                    
                end
                
                for i=1:N_enabled_thruster
                    if thruster_data(i).type == 3
                        [phi_0(i),phiPlus0(i),inside_fz(i)] = angle_point_to(...
                            alloc_out(i).phi, ...
                            thruster_data(i).phi_reserve, ...
                            thruster_data(i).dphi_max, ...
                            dt, ...
                            thruster_data(i).phi_min, ...
                            thruster_data(i).phi_max);
                        on_border(i) = inside_fz(i);
                        if inside_fz(i) == 0
                            phi_guider_old(i) = alloc_out(i).phi;
                            
                        end
                        if inside_fz(i) == 0%phi is outside forbidden zone
                            [thruster_data(i).phi_,thruster_data(i).phiPlus]= ...
                                angleMaxMin(thruster_data(i).phi_min, ...
                                thruster_data(i).phi_max, ...
                                thruster_data(i).phi_reserve, ...
                                phi_0(i), ...
                                phiPlus0(i));
                            thruster_data(i).T_ = max(thruster_data(i).T_reserve - ...
                                dt * thruster_data(i).dTmax, thruster_data(i).Tmin);
                            thruster_data(i).Tplus = min(thruster_data(i).T_reserve + ...
                                dt * thruster_data(i).dTmax, thruster_data(i).Tmax);
                        else
                            thruster_data(i).phi_ = phi_0(i);
                            thruster_data(i).phiPlus = phiPlus0(i);
                            thruster_data(i).T_ = thruster_data(i).Tmin;
                            thruster_data(i).Tplus = thruster_data(i).Tmin + 10;
                            
                        end
                        thruster_data(i).phi = thruster_data(i).phi_ ...
                            + angle_dist(thruster_data(i).phi_,thruster_data(i).phiPlus) / 2;
                        
                        thruster_data(i).T = thruster_data(i).T_reserve;
                        
                        
%                         if abs(thruster_data(i).phi_max - thruster_data(i).phi_min)> FeasibilityTol%means there is forbidden zone
%                         
%                            if min([angle_dist(alloc_out(i).phi,thruster_data(i).phi_min)
%                                    angle_dist(thruster_data(i).phi_min,alloc_out(i).phi)
%                                    angle_dist(alloc_out(i).phi,thruster_data(i).phi_max)
%                                    angle_dist(thruster_data(i).phi_max,alloc_out(i).phi)]) ...
%                                    <= FeasibilityTol
%                                on_border(i) = 1;
%                            else
%                                on_border(i) = 0;
%                            end
%                            
%                         
%                         end
                        
                        
                    else
                        thruster_data(i).phi = thruster_data(i).phi_reserve;
                        thruster_data(i).T = thruster_data(i).T_reserve;
                    end
                    
                end
            end
            
            
        else  %if step == 2
            for i=1:N_enabled_thruster
                thruster_data(i).phi_reserve = alloc_out(i).phi;%previous phi
                thruster_data(i).T_reserve = alloc_out(i).T;%previous T
                thruster_data(i).phi = alloc_out(i).phi;
                thruster_data(i).T = alloc_out(i).T;
            end
            run = false;%break the while loop
        end
        
    end
end

% status
% chck_sol = [[thruster_data(:).phi_min]' [thruster_data(:).phi_]' [alloc_out.phi]' [thruster_data.phiPlus]' [thruster_data(:).phi_max]'];
% check_angle(chck_sol(:,1),chck_sol(:,2),chck_sol(:,3),chck_sol(:,4),chck_sol(:,5))
% solution
end

function [T,new,in_fb] = baseThrust(T,phi,T_base,phi_base,phi_min,phi_max)
tol = 1e-6;
x = T * cos(phi) + T_base * cos(phi_base);
y = T * sin(phi) + T_base * sin(phi_base);


T = norm([x y]);
new = atan2(y,x);
in_fb = 0;
if abs(phi_max - phi_min)> tol%means there is forbidden zone
if angle_dist(phi_min,new) > angle_dist(phi_min,phi_max) %inside fbzone
    in_fb = 1;
end
end
end

function [phi_0,phiPlus0,inside_fz] = angle_point_to(s,phi,dphi_max,dt,phi_min,phi_max)
tol = 1e-6;
%s: setpoint from global solution
inside_fz = 0;%phi, phi_0 and phiPlus0 is outside forbidden zone
if min(angle_dist(phi,s),angle_dist(s,phi)) > dphi_max * dt
    move = dphi_max * dt / 2;
else
    move = 0;
end

% find out if phi, phi_0 or phiPlus0 is inside fb zone
if angle_dist(phi,s) < pi
    phi_0 = phi + move;
    phiPlus0 = phi + dphi_max * dt;
    
else
    phiPlus0 = phi - move;
    phi_0 = phi - dphi_max * dt;
end
if abs(phi_max - phi_min)> tol%means there is forbidden zone
if angle_dist(phi_min,phi) >= angle_dist(phi_min,phi_max)% ...
%         || angle_dist(phi_min,phi_0) >= angle_dist(phi_min,phi_max) ...
%         || angle_dist(phi_min,phiPlus0) >= angle_dist(phi_min,phi_max)
%       commented at v4.10
    inside_fz = 1;%phi, phi_0 or phiPlus0 is inside forbidden zone
end
end
% % % % % % % % %
if abs(phi_max - phi_min)> tol%means there is forbidden zone
    if inside_fz == 0 %phi is outside forbidden zone
        if angle_dist(phi,s) < pi
            phi_0 = phi + move;
            phiPlus0 = phi + dphi_max * dt;
        else
            phiPlus0 = phi - move;
            phi_0 = phi - dphi_max * dt;
        end
    else%phi is inside forbidden zone
        if angle_dist(phi_min,s) <= angle_dist(phi_min,phi_max)%s is outside forbidden zone
            if angle_dist(phi,s) < pi
                phi_0 = phi + move;
                phiPlus0 = phi + dphi_max * dt;
            else
                phiPlus0 = phi - move;
                phi_0 = phi - dphi_max * dt;
            end
        else% s is inside forbidden zone
            if angle_dist(phi_max,phi) < angle_dist(phi,phi_min) %phi is closer to phi_max
                phiPlus0 = phi - dphi_max * dt / 2;
                phi_0 = phi - dphi_max * dt;
            else%phi is closer to phi_min
                phi_0 = phi + dphi_max * dt / 2;
                phiPlus0 = phi + dphi_max * dt;
            end
        end
    end
else%there is no forbidden zone
    if angle_dist(phi,s) < pi
        phi_0 = phi + move;
        phiPlus0 = phi + dphi_max * dt;
    else
        phiPlus0 = phi - move;
        phi_0 = phi - dphi_max * dt;
    end
end
end
