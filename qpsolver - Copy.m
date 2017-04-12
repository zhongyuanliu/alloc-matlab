    
function [solution,status]=qpsolver(thruster_data,alloc_out,T_r,N_enabled_thruster)
coder.inline('never');

Ntunnel=0;
Nazi=0;
dt = thruster_data(1).dt;
c2 = 0.001;
% ---------------precalculation------------------
for j=1:N_enabled_thruster
    switch thruster_data(j).type
        case 4%tunnel
            
            thruster_data(j).Tplus = min(thruster_data(j).T + ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmax);
            thruster_data(j).T_ = max(thruster_data(j).T - ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmin);  
            
            Ntunnel = Ntunnel + 1;
        case 3%azi
            thruster_data(j).Tplus = min(thruster_data(j).T + ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmax);
            thruster_data(j).T_ = max(thruster_data(j).T - ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmin);
            thruster_data(j).phiPlus = min(thruster_data(j).phi + ...
                dt * thruster_data(j).dphi_max, thruster_data(j).phi_max);
            thruster_data(j).phi_ = max(thruster_data(j).phi - ...
                dt * thruster_data(j).dphi_max, thruster_data(j).phi_min);
            Nazi = Nazi + 1;
    end
end
Nvar = 2 * Ntunnel + 2 * Nazi + 3;%number of variables in quad
N_A = Ntunnel * 2 + Nazi * 5;
% ------------------------------------------------------
H = eye(Nvar);
f= zeros(Nvar,1);
N_Aeq = 3 + Ntunnel;
Aeq= zeros(N_Aeq,Nvar);
beq=zeros(N_Aeq,1);
beq(1) = T_r.Tx;
beq(2) = T_r.Ty;
beq(3) = T_r.Tm;
A=zeros(N_A,Nvar);
b= zeros(N_A,1);
% ------------------------------------------------------
i_H = 1;
numOfTunnel=0;
numOfAzi=0;
pos_A = 1;
pos_Aeq = 4;
for j=1:N_enabled_thruster
    switch thruster_data(j).type
        case 4  %tunnel
            H(i_H,i_H) = thruster_data(j).weight;
            H(i_H + 1,i_H + 1) = thruster_data(j).weight;
            i_H = i_H + 2;
            
            for i = 1:3%R1 R2 R3 
                switch i
                    case 1
                        Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + 1) = 1;
                        Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + 2) = 0;
                    case 2
                        Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + 1) = 0;
                        Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + 2) = 1;
                    case 3
                        Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + 1) = ...
                            - thruster_data(j).y + thruster_data(j).y0;
                        Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + 2) = ...
                            thruster_data(j).x - thruster_data(j).x0;
                end

            end
                
            Aeq(pos_Aeq,numOfAzi * 2 + numOfTunnel * 2 + 1) = sin(thruster_data(j).phi);
            Aeq(pos_Aeq,numOfAzi * 2 + numOfTunnel * 2 + 2) = -cos(thruster_data(j).phi);
            beq(pos_Aeq) = 0;                 
           
            
            A(pos_A,numOfAzi * 2 + numOfTunnel * 2 + 1) = -cos(thruster_data(j).phi);
            A(pos_A,numOfAzi * 2 + numOfTunnel * 2 + 2) = -sin(thruster_data(j).phi);
            b(pos_A) = -thruster_data(j).Tplus;
            A(pos_A + 1,numOfAzi * 2 + numOfTunnel * 2 + 1) = cos(thruster_data(j).phi);
            A(pos_A + 1,numOfAzi * 2 + numOfTunnel * 2 + 2) = sin(thruster_data(j).phi);
            b(pos_A + 1) = thruster_data(j).T_;
            
            pos_A = pos_A + 2;
            pos_Aeq = pos_Aeq + 1;
            numOfTunnel = numOfTunnel + 1;
        case 3  %azi
            H(i_H,i_H) = thruster_data(j).weight;
            H(i_H + 1,i_H + 1) = thruster_data(j).weight;
            i_H = i_H + 2;
            for ixy =1:2 % xthurst and ythrust
                for i = 1:3 %R1 R2 R3
                    if ixy == 1
                        switch i
                            case 1
                                Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + ixy) = 1;
                            case 2
                                Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + ixy) = 0;
                            case 3
                                Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + ixy) = ...
                                    - thruster_data(j).y + thruster_data(j).y0;

                        end
                        
                    elseif ixy ==2
                        switch i
                            case 1
                                Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + ixy) = 0;
                            case 2
                                Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + ixy) = 1;
                            case 3
                                Aeq(i,numOfAzi * 2 + numOfTunnel * 2 + ixy) = ...
                                    thruster_data(j).x - thruster_data(j).x0;

                        end                        
                    end
                end
                
                if ixy == 1
                    A(pos_A,numOfAzi * 2 + numOfTunnel * 2 + ixy) = -sin(thruster_data(j).phi_);%(3.15a)
                    b(pos_A) = -c2 * abs(cos(thruster_data(j).phi_)) * ...
                        (thruster_data(j).Tplus - thruster_data(j).T_);
                    
                    A(pos_A + 1,numOfAzi * 2 + numOfTunnel * 2 + ixy) = sin(thruster_data(j).phiPlus);%(3.15b)
                    b(pos_A + 1) = -c2 * abs(cos(thruster_data(j).phiPlus)) * ...
                        (thruster_data(j).Tplus - thruster_data(j).T_);
                    
                    A(pos_A + 2,numOfAzi * 2 + numOfTunnel * 2 + ixy) = cos(thruster_data(j).phi);%(3.16a)
                    b(pos_A + 2) = thruster_data(j).T_;
                    
                    A(pos_A + 3,numOfAzi * 2 + numOfTunnel * 2 + ixy) = -cos(thruster_data(j).phi_);%(3.16b)
                    b(pos_A + 3) = -cos(thruster_data(j).phi_ - thruster_data(j).phi) * ...
                        thruster_data(j).Tplus;
                     
                    A(pos_A + 4,numOfAzi * 2 + numOfTunnel * 2 + ixy) = -cos(thruster_data(j).phiPlus);%(3.16c)
                    b(pos_A + 4) = -cos(thruster_data(j).phiPlus - thruster_data(j).phi) * ...
                        thruster_data(j).Tplus;
                    
                elseif ixy == 2
                    A(pos_A,numOfAzi * 2 + numOfTunnel * 2 + ixy) = cos(thruster_data(j).phi_);%(3.15a)
                                  
                    A(pos_A + 1,numOfAzi * 2 + numOfTunnel * 2 + ixy) = -cos(thruster_data(j).phiPlus);%(3.15b)
                   
                    A(pos_A + 2,numOfAzi * 2 + numOfTunnel * 2 + ixy) = sin(thruster_data(j).phi);%(3.16a)
                    
                    A(pos_A + 3,numOfAzi * 2 + numOfTunnel * 2 + ixy) = -sin(thruster_data(j).phi_);%(3.16b)
                    
                    A(pos_A + 4,numOfAzi * 2 + numOfTunnel * 2 + ixy) = -sin(thruster_data(j).phiPlus);%(3.16c)
                   
                end
                
            end
            pos_A = pos_A + 5;

            

            numOfAzi = numOfAzi + 1;
    end
end

Aeq(1,numOfTunnel * 2 + 2 * numOfAzi + 1) = -1;
Aeq(2,numOfTunnel * 2 + 2 * numOfAzi + 2) = -1;
Aeq(3,numOfTunnel * 2 + 2 * numOfAzi + 3) = -1;

for i = 1:3
H(i_H + i - 1,i_H + i - 1) = thruster_data(1).weight_s(i);%!!!thruster_data(1) must be available
end

MaxIter=100;
FeasibilityTol=1e-5;

[L,p] = Chol_fc(H);
% before chol

if p==0
Linv = inv(L);
opt = mpcqpsolverOptions;
opt.MaxIter = MaxIter;
opt.IntegrityChecks = false;
opt.FeasibilityTol = FeasibilityTol;
iA0 = false(size(b));
[solution,status,iA,lambda] = Mpcqpsolver(Linv,f,A,b,Aeq,beq,iA0,opt);
% x = quadprog(H,f,A,b,Aeq,beq);
else
  solution = zeros(Nvar,1);
  status = -2;
end

% output
numOfTunnel = 0;
numOfAzi = 0;
for i=1:N_enabled_thruster
    alloc_out(i).label = thruster_data(i).label;
    alloc_out(i).enable = thruster_data(i).enable;
    
    switch thruster_data(i).type
        case 4
            alloc_out(i).Tx = 0;
            alloc_out(i).Ty = solution(numOfTunnel * 2 + numOfAzi * 2 + 1);
            alloc_out(i).T = alloc_out(i).Ty;
            alloc_out(i).Tm = alloc_out(i).Ty * (thruster_data(i).x ...
                - thruster_data(i). x0) + alloc_out(i).Tx ...
                * (thruster_data(i).y0 - thruster_data(i).y);
            
            numOfTunnel = numOfTunnel + 1;
        case 3
            alloc_out(i).Tx = solution(numOfTunnel * 2 + numOfAzi * 2 + 1);
            alloc_out(i).Ty = solution(numOfTunnel * 2 + numOfAzi * 2 + 2);
            alloc_out(i).T = sqrt(alloc_out(i).Tx^2 + alloc_out(i).Ty^2);
            alloc_out(i).phi = atan2(solution(numOfTunnel * 2 + numOfAzi * 2 + 2), ...
                solution(numOfTunnel * 2 + numOfAzi * 2 + 1));
            alloc_out(i).Tm = alloc_out(i).Ty * (thruster_data(i).x ...
                - thruster_data(i). x0) + alloc_out(i).Tx ...
                * (thruster_data(i).y0 - thruster_data(i).y);
           numOfAzi = numOfAzi + 1; 
    end
            
            
end

  