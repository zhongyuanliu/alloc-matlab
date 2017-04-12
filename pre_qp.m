function [thruster_data,Nvar,Ho,fo,Ao,bo,Aeqo,beqo,x0] = pre_qp(init,thruster_data,T_r,N_enabled_thruster,Ce,NA_fpp,rudder_table,Rangle_T,Rangle_Fangle,no_azi_angle_constr,method)
%  1.04 thruster_data(i).phi is replaced by thruster_data(i).phi_min for tunnel thruster
%  1.06 guarantee initial x feasible
%  4.06 u is scaled. c_scale * Tmax^2 == 1 |13-03-2017 by Liu
%% ---------------precalculation------------------
persistent A_azi NA_azi A_fpp b_fpp;
tol = 1e-6;%!!!avoid numerical problem
Ntunnel=0;
Nazi=0;
Nfpp=0;
dt = thruster_data(1).dt;
c2 = 0.001*0;%!!!!!!!!!!!!!!!!!!!!!!!!!!

error_linearization = 0.05;
if isempty(A_azi) || isempty(NA_azi) || init == true
    [A_azi,NA_azi] = N_sided_linear(error_linearization);
end

if isempty(A_fpp) || isempty(b_fpp) || init == true
    [A_fpp,b_fpp,~] = N_sided_rudder(1e5 * rudder_table(:,2:3));%%!!what if there is no fpp!!!!!!!!!!
end
% NA_azi = 0;%not used at the moment!!!!!!!!!!!!!!!!

% init = true;
x0=zeros(N_enabled_thruster * 2 + 3,1);

for j=1:N_enabled_thruster
%     check if phi meets constraint, if check_phi == 0, phi is not between
%     max and min, then it should be changed. VERY IMPORTANT
%     check_phi = check_angle(thruster_data(j).phi_min,thruster_data(j).phi_min, ...
%         thruster_data(j).phi,thruster_data(j).phi_max,thruster_data(j).phi_max);
%     
%     if check_phi == 0 % guarantee initial phi feasible
%         thruster_data(j).phi = ...
%             thruster_data(j).phi_min + ...
%             angle_dist(thruster_data(j).phi_min, thruster_data(j).phi_max) / 2;
% 
%         thruster_data(j).phi = Pi_toPi(thruster_data(j).phi);
%     end
    
    if thruster_data(j).type == 4%added, avoid tunnel infeasible angle input, v1.07
        if abs(thruster_data(j).phi_min - thruster_data(j).phi) > tol
            thruster_data(j).phi = thruster_data(j).phi_min;
        end
        
    end
    
    if abs(thruster_data(j).Tmax - thruster_data(j).Tmin) < tol%avoid Tmax == Tmin
        thruster_data(j).Tmax = 1e7;%big enough
        if thruster_data(j).type == 4
            thruster_data(j).Tmin = -1e7;
        else
            thruster_data(j).Tmin = 0;
        end

    end
    
    if thruster_data(j).T > thruster_data(j).Tmax ...
            || thruster_data(j).T < thruster_data(j).Tmin% guarantee initial T feasible 
        if thruster_data(j).type == 4
            thruster_data(j).T = 0;
        else
            thruster_data(j).T = thruster_data(j).Tmin;
        end
    end
    
    if thruster_data(j).T_reserve > thruster_data(j).Tmax ...
            || thruster_data(j).T_reserve < thruster_data(j).Tmin% guarantee initial T feasible 
        if thruster_data(j).type == 4
            thruster_data(j).T_reserve = 0;
        else
            thruster_data(j).T_reserve = thruster_data(j).Tmin;
        end
    end    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    
    x0(2 * j - 1)=thruster_data(j).T*cos(thruster_data(j).phi);
    x0(2 * j)=thruster_data(j).T*sin(thruster_data(j).phi);
    
    switch thruster_data(j).type
        case 4%tunnel
            
            thruster_data(j).Tplus = min(thruster_data(j).T + ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmax);
            
            thruster_data(j).T_ = max(thruster_data(j).T - ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmin);
            
            Ntunnel = Ntunnel + 1;
        case 3%azi
%             thruster_data(j).Tplus = min(thruster_data(j).T + ...
%                 dt * thruster_data(j).dTmax, thruster_data(j).Tmax);
%             
%             thruster_data(j).T_ = max(thruster_data(j).T - ...
%                 dt * thruster_data(j).dTmax, thruster_data(j).Tmin);
%             -----------------------------
%             if thruster_data(j).constr == 2 
%                 thruster_data(j).phi_ = thruster_data(j).phi;%!!!!!!add -tol or not
%                 thruster_data(j).phiPlus = thruster_data(j).phi;%!!!!!add +tol or not
%             elseif thruster_data(j).constr == 3 
%                 
%             end
%             thruster_data(j).phi_ = thruster_data(j).phi - dt * thruster_data(j).dphi_max;
%                 thruster_data(j).phiPlus = thruster_data(j).phi + dt * thruster_data(j).dphi_max;
%             -----------------------------
            %             thruster_data(j).phiPlus = min(thruster_data(j).phi + ...
            %                 dt * thruster_data(j).dphi_max, thruster_data(j).phi_max);
            %
            %             thruster_data(j).phi_ = max(thruster_data(j).phi - ...
            %                 dt * thruster_data(j).dphi_max, thruster_data(j).phi_min);
            
            Nazi = Nazi + 1;
        case 2%fpp
            thruster_data(j).Tplus = min((thruster_data(j).T + ...
                dt * thruster_data(j).dTmax) , ...
                thruster_data(j).Tmax) * polyval(Rangle_T,abs(thruster_data(j).phi));
            
            thruster_data(j).T_ = max(thruster_data(j).T - ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmin) ...
                * polyval(Rangle_T,abs(thruster_data(j).phi));
            
            
            [temp1,temp2]= ...
                angleMaxMin(thruster_data(j).phi_min, ...
                thruster_data(j).phi_max, ...
                thruster_data(j).phi, ...
                thruster_data(j).phi - dt * thruster_data(j).dphi_max, ...
                thruster_data(j).phi + dt * thruster_data(j).dphi_max);
            thruster_data(j).phiPlus = sign(temp2) * polyval(Rangle_Fangle,abs(temp2));
            thruster_data(j).phi_ = sign(temp1) * polyval(Rangle_Fangle,abs(temp1));
            %             temp = (min(thruster_data(j).phi + ... %force angle
            %                 dt * thruster_data(j).dphi_max, thruster_data(j).phi_max));
            %             thruster_data(j).phiPlus = sign(temp) * polyval(Rangle_Fangle,abs(temp));
            
            %             temp = (max(thruster_data(j).phi - ... %force angle
            %                 dt * thruster_data(j).dphi_max, thruster_data(j).phi_min));
            %             thruster_data(j).phi_ = sign(temp) * polyval(Rangle_Fangle,abs(temp));
            
            Nfpp = Nfpp + 1;
    end
end

    
%!!!!!to be improved

%     for j=1:N_enabled_thruster%!!!!!the above Tplus is overwitten
        
%         thruster_data(j).Tplus = thruster_data(j).Tmax;
%         thruster_data(j).T_ = thruster_data(j).Tmin;
%         thruster_data(j).phi_ = thruster_data(j).phi_min;%not correct. if phi_max - phi_min >pi, constraint falls
%         thruster_data(j).phiPlus = thruster_data(j).phi_max;
        
%     end

% chck_sol = [[thruster_data(:).phi_min]' [thruster_data(:).phi_]' [thruster_data.phi]' [thruster_data.phiPlus]' [thruster_data(:).phi_max]'];
% check_angle(chck_sol(:,1),chck_sol(:,2),chck_sol(:,3),chck_sol(:,4),chck_sol(:,5))

Nvar = 2 * Ntunnel + 2 * Nazi + 2 * Nfpp + 3;%number of variables in quad
% if no_azi_angle_constr == 0 % there is angle constraint
%     N_A = Ntunnel * 2 + Nazi * (5 + NA_azi) + Nfpp * (5 + NA_fpp);
% else
%     N_A = Ntunnel * 2 + Nazi *  NA_azi + Nfpp * (5 + NA_fpp);
% end

N_A = 0;%added v1.08
for j=1:N_enabled_thruster %added v1.08
    switch thruster_data(j).type
        case 4
            N_A = N_A + 2;
        case 3

            if thruster_data(j).constr == 2 || thruster_data(j).constr == 3
                N_A = N_A + (5 + NA_azi);
            else 
                N_A = N_A + NA_azi;
            end
        case 2
            N_A = N_A + (5 + NA_fpp);
    end
    
end
%% ------------------------------------------------------
persistent H f N_Aeq Aeq beq;
if isempty(H)||isempty(f)||isempty(N_Aeq)||isempty(Aeq)||isempty(beq)||init == true
    H = eye(Nvar);
    f= zeros(Nvar,1);
    N_Aeq = 3 + Ntunnel;
    Aeq= zeros(N_Aeq,Nvar);
    beq=zeros(N_Aeq,1);


end
    A=zeros(N_A,Nvar);
    b= zeros(N_A,1);
    beq(1) = T_r.Tx;
    beq(2) = T_r.Ty;
    beq(3) = T_r.Tm;
% ------------------------------------------------------
i_H = 1;
numOfTunnel=0;
numOfAzi=0;
numOfFpp=0;
pos_A = 1;
pos_Aeq = 4;

eff_all = ones(1,N_enabled_thruster);
for j=1:N_enabled_thruster
    if thruster_data(j).type == 3
        eff_all(j) = azi_eff(Ce(j).value,thruster_data(j).phi);
    end    
end
[eff,pos_eff] = min(eff_all);
eff_all = ones(1,N_enabled_thruster);
eff_all(pos_eff) = eff;

for j=1:N_enabled_thruster
    switch thruster_data(j).type
        case 4  %tunnel
            temp = numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2;
            if init == true
            H(i_H,i_H) = thruster_data(j).weight(1) * thruster_data(j).Tmax^4;
            H(i_H + 1,i_H + 1) = thruster_data(j).weight(2) * thruster_data(j).Tmax^4;
            
            

            
                Aeq(1:3,temp + 1:temp + 2)...
                    = [1 0;0 1; - thruster_data(j).y + thruster_data(j).y0 thruster_data(j).x - thruster_data(j).x0] * thruster_data(j).Tmax^2;%R1 R2 R3
                
                
%             Aeq(pos_Aeq,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = sin(thruster_data(j).phi_min);
%             Aeq(pos_Aeq,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2) = -cos(thruster_data(j).phi_min);
            Aeq(pos_Aeq,temp + 1:temp + 2) = [sin(thruster_data(j).phi_min) -cos(thruster_data(j).phi_min)] * thruster_data(j).Tmax^2;
                
            
            beq(pos_Aeq) = 0;
            
            
            end
            %             A(pos_A,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = -cos(thruster_data(j).phi_min);
%             A(pos_A,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2) = -sin(thruster_data(j).phi_min);
            A(pos_A,temp + 1:temp + 2) = [-cos(thruster_data(j).phi_min) -sin(thruster_data(j).phi_min)]* thruster_data(j).Tmax^2;
            
%             b(pos_A) = -thruster_data(j).Tplus;
%             A(pos_A + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = cos(thruster_data(j).phi_min);
%             A(pos_A + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2) = sin(thruster_data(j).phi_min);
            A(pos_A + 1,temp + 1:temp + 2) = [cos(thruster_data(j).phi_min) sin(thruster_data(j).phi_min)]* thruster_data(j).Tmax^2;

%             b(pos_A + 1) = thruster_data(j).T_;
            b(pos_A:pos_A + 1) = [-thruster_data(j).Tplus; thruster_data(j).T_];
            %             check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 1, ...
            %                 numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ...
            %                 b(pos_A:pos_A + 1));
            i_H = i_H + 2;
            pos_A = pos_A + 2;
            pos_Aeq = pos_Aeq + 1;
            numOfTunnel = numOfTunnel + 1;
        case 3  %azi
            if init == true
            H(i_H,i_H) = thruster_data(j).weight(1) * thruster_data(j).Tmax^4;
            H(i_H + 1,i_H + 1) = thruster_data(j).weight(2) * thruster_data(j).Tmax^4;
            end
            
            temp = numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2;
            Aeq(1:3,temp + 1:temp + 2) = [eff_all(j) 0;0 eff_all(j); - thruster_data(j).y + thruster_data(j).y0 thruster_data(j).x - thruster_data(j).x0] * thruster_data(j).Tmax^2;
            if thruster_data(j).constr == 2 || thruster_data(j).constr == 3 % there is constraint
                A(pos_A:pos_A + 4,temp + 1 : temp + 2) = [-sin(thruster_data(j).phi_) cos(thruster_data(j).phi_)
                    sin(thruster_data(j).phiPlus) -cos(thruster_data(j).phiPlus)
                    cos(thruster_data(j).phi) sin(thruster_data(j).phi)
                    -cos(thruster_data(j).phi_) -sin(thruster_data(j).phi_)
                    -cos(thruster_data(j).phiPlus) -sin(thruster_data(j).phiPlus)] * thruster_data(j).Tmax^2;
                b(pos_A : pos_A + 4) = [-c2 * abs(cos(thruster_data(j).phi_)) * ...
                    (thruster_data(j).Tplus - thruster_data(j).T_)
                    -c2 * abs(cos(thruster_data(j).phiPlus)) * ...
                    (thruster_data(j).Tplus - thruster_data(j).T_)
                    thruster_data(j).T_
                    -cos(thruster_data(j).phi_ - thruster_data(j).phi) * ...
                    thruster_data(j).Tplus
                    -cos(thruster_data(j).phiPlus - thruster_data(j).phi) * ...
                    thruster_data(j).Tplus];
%                 if init == true
                    A(pos_A + 4 + 1:pos_A + 4 + NA_azi,temp + 1:temp + 2) = A_azi * thruster_data(j).Tmax^2;
                    b(pos_A + 4 + 1:pos_A + 4 + NA_azi) = - thruster_data(j).Tmax * cos(pi / NA_azi);
%                 end
            else
%                 if init == true
                    A(pos_A - 1 + 1 : pos_A - 1 + NA_azi,temp + 1:temp + 2) = A_azi * thruster_data(j).Tmax^2;
                    
                    b(pos_A - 1 + 1 : pos_A - 1 + NA_azi) = - thruster_data(j).Tmax * cos(pi / NA_azi);
%                 end
            end    
                
            %             check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 4, ...
            %                 numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ...
            %                 b(pos_A:pos_A + 4));%azi
            %
            %                     check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A + 4 + 1:pos_A + 4 + i_Azi, ...
            %                         numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ...
            %                         b(pos_A + 4 + 1:pos_A + 4 + i_Azi));
            %         check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 4 + i_Azi, ...
            %             numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ...
            %             b(pos_A:pos_A + 4 + i_Azi));
            if thruster_data(j).constr == 2 || thruster_data(j).constr == 3 % there is constraint 
                pos_A = pos_A + 5 + NA_azi;
            else
                pos_A = pos_A + NA_azi;
            end
            
            
            i_H = i_H + 2;
            numOfAzi = numOfAzi + 1;
        case 2 %fpp
            [A_fpp,b_fpp,~] = N_sided_rudder(thruster_data(j).Tmax * rudder_table(:,2:3));%x,y
            if init == true
            H(i_H,i_H) = thruster_data(j).weight(1) * thruster_data(j).Tmax^4;
            H(i_H + 1,i_H + 1) = thruster_data(j).weight(2) * thruster_data(j).Tmax^4;
            end
            temp = numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2;
            if init == true
            Aeq(1:3,temp + 1:temp + 2) = [1 0;0 1; - thruster_data(j).y + thruster_data(j).y0 thruster_data(j).x - thruster_data(j).x0] * thruster_data(j).Tmax^2;
            end
            A(pos_A:pos_A + 4,temp + 1 : temp + 2) = [-sin(thruster_data(j).phi_) cos(thruster_data(j).phi_)
                    sin(thruster_data(j).phiPlus) -cos(thruster_data(j).phiPlus)
                    cos(thruster_data(j).phi) sin(thruster_data(j).phi)
                    -cos(thruster_data(j).phi_) -sin(thruster_data(j).phi_)
                    -cos(thruster_data(j).phiPlus) -sin(thruster_data(j).phiPlus)] * thruster_data(j).Tmax^2;
                b(pos_A : pos_A + 4) = [-c2 * abs(cos(thruster_data(j).phi_)) * ...
                    (thruster_data(j).Tplus - thruster_data(j).T_)
                    -c2 * abs(cos(thruster_data(j).phiPlus)) * ...
                    (thruster_data(j).Tplus - thruster_data(j).T_)
                    thruster_data(j).T_
                    -cos(thruster_data(j).phi_ - thruster_data(j).phi) * ...
                    thruster_data(j).Tplus
                    -cos(thruster_data(j).phiPlus - thruster_data(j).phi) * ...
                    thruster_data(j).Tplus];
                if init == true
                 A(pos_A + 4 + 1:pos_A + 4 + NA_fpp,temp + 1:temp + 2) = A_fpp * thruster_data(j).Tmax^2;
                    b(pos_A + 4 + 1:pos_A + 4 + NA_fpp) = b_fpp;
                end
            
            %         check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 4, ...
            %             numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ...
            %             b(pos_A:pos_A + 4));%fpp
            %
            %         check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A + 4 + 1:pos_A + 4 + i_Fpp, ...
            %             numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ...
            %             b(pos_A + 4 + 1:pos_A + 4 + i_Fpp));
            %             check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 4 + i_Fpp, ...
            %                 numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ...
            %                 b(pos_A:pos_A + 4 + i_Fpp));
            i_H = i_H + 2;
            pos_A = pos_A + 5 + NA_fpp;
            
            numOfFpp = numOfFpp + 1;
    end
end
%!!-----------------add bound for s1 s2 s3---does not help------------
% A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = 1;
% A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = 1;
% A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = 1;
% A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = -1;
% A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = -1;
% A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = -1;
% b(end + 1) = -1000000;
% b(end + 1) = -1000000;
% b(end + 1) = -1000000;
% b(end + 1) = 1000000;
% b(end + 1) = 1000000;
% b(end + 1) = 1000000;
%-------------------------------------------------
if init == true
temp = numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2;
    Aeq(1,temp + 1) = -1;%s1
    Aeq(2,temp + 2) = -1;%s2
    Aeq(3,temp + 3) = -1;%s3
for i = 1:3
    H(i_H + i - 1,i_H + i - 1) = thruster_data(1).weight_s(i);%!!!thruster_data(1) must be available
end
end
Ho = H;
fo = f;
Ao = A;
bo = b;
Aeqo = Aeq;
beqo = beq;
end

function [A,N] = N_sided_linear(error)
%%
N = floor(pi / (acos(1 - error))) + 1;
%     r = - R * cos(pi / N);%change <= to >=
phi = zeros(N,1);
A = zeros(N,2);
for i = 0 : (N - 1)
    phi(i + 1) = (2 * i + 1) * pi / N;
    A(i + 1,1) = - cos(phi(i + 1));%change <= to >=
    A(i + 1,2) = - sin(phi(i + 1));%change <= to >=
end

end

function [A,b,N] = N_sided_rudder(points)%x,y
%%
% points =[0.980,0;0.965,0.0800;0.956,0.163;0.920,0.200;0.895,0.287;0.840,0.300
%     0.811,0.341;0.770,0.345;0.743,0.347;0.720,0.300];
temp = flipud(points(2:end,:));
temp(:,2) = -temp(:,2);
points = [0 0;temp;points;0 0];

N = size(points,1) - 1;
b = zeros(N,1);
A = zeros(N,2);
%     points(end + 1,:) = points(1,:);
for i = 1 : N
    A(i,1) = -points(i + 1,2) + points(i,2);
    A(i,2) = -points(i,1) + points(i + 1,1);
    b(i) = -points(i,1) * points(i + 1,2) + points(i + 1,1) * points(i,2);
    
end

% check_rudder_constr(points,A,b);
end

function check = check_rudder_constr(points,A,b)
%%
check = 0;
[max_dat]=max(points,[],1);
[min_dat]=min(points,[],1);
interv = (max_dat - min_dat)/30;
[x, y]=meshgrid(min_dat(1):interv:max_dat(1),min_dat(2):interv:max_dat(2));
M=size(x,1);
NN=size(x,2);
dat_sav = [];

plot(points(:,2),points(:,1));
axis equal
hold on
for i = 1:M
    for j=1:NN
        if A * [x(i,j) y(i,j)]' >= b
            dat_sav(end+1,:) = [x(i,j) y(i,j)];
            plot(y(i,j),x(i,j),'*');
            check = 1;
        end
    end
end
hold off
end

function y = check_angle(phi_min,phi_,phi,phiPlus,phi_max)
tol = 1e-8;%!!!avoid numerical problem
if abs(phi_min - phi_max) < tol%phi_min == phi_max means no forbidden zone
    k = true;
else
k = (angle_dist(phi_min,phi_)<=angle_dist(phi_min,phi) + tol) ...
    &(angle_dist(phi_min,phi)<=angle_dist(phi_min,phiPlus) + tol) ...
    & (angle_dist(phi,phiPlus)<=angle_dist(phi,phi_max) + tol);
end
if prod(k) == 0
    y = 0;
else
    y = 1;
end

end



