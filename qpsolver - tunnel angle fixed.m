    
function [solution,status]=qpsolver(thruster_data,alloc_out,T_r,N_enabled_thruster)
coder.inline('never');

% alloc_out = struct ('label',{1,2,3,4,5,6,7,8}, ...
%                     'enable',{1,1,1,1,1,1,1,1}, ...
%                     'TSP',0, ...
%                     'ASP',0, ...
%                     'Tx',0, ...
%                     'Ty',0, ...
%                     'T',0, ...
%                     'phi',0, ...
%                     'Tm',0);
% thruster_data = struct('label',{1,2,3,4,5,6,7,8}, ...
%                         'type',{3,3,3,3,3,3,3,3}, ...
%                         'enable',{1,1,1,1,1,1,1,1}, ...
%                         'dt',0.5, ...
%                         'x0',0, ...
%                         'y0',0, ...
%                         'x',{10,-30,-30,21,23,18,5,5}, ...
%                         'y',{4,-3.5,3.5,0,0,0,4,-4}, ...
%                         'T',{5435,2305,351,8955,3211,231,25513,5866}, ...
%                         'phi',{1.02,0.456,3,1.432,2.8,0.2,0.08,3.12}, ...
%                         'Tmax',71500, ...
%                         'Tmin',1000, ...
%                         'dTmax',7150, ...
%                         'dphi_max',0.5, ...
%                         'Tplus',0, ...
%                         'T_',0, ...
%                         'phiPlus',0, ...
%                         'phi_',0, ...
%                         'phi_min',{-1.04,-1.04,2.09,0.7,0.1,-2.4,-2.4,0.01}, ...
%                         'phi_max',{-2.09,-2.09,1.04,-0.7,6,2.4,2.4,6}, ...
%                         'weight',1,...
%                         'weight_s',[1000000 1000000 100000000]');  
thruster_data
numOfTunnel=0;
numOfAzi=0;
dt = thruster_data(1).dt;
c2 = 0.001;
%precalculation
for j=1:N_enabled_thruster
    switch thruster_data(j).type
        case 4%tunnel
            
            thruster_data(j).Tplus = min(thruster_data(j).T + ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmax);
            thruster_data(j).T_ = max(thruster_data(j).T - ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmin);  
            
            numOfTunnel = numOfTunnel + 1;
        case 3%azi
            thruster_data(j).Tplus = min(thruster_data(j).T + ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmax);
            thruster_data(j).T_ = max(thruster_data(j).T - ...
                dt * thruster_data(j).dTmax, thruster_data(j).Tmin);
            thruster_data(j).phiPlus = min(thruster_data(j).phi + ...
                dt * thruster_data(j).dphi_max, thruster_data(j).phi_max);
            thruster_data(j).phi_ = max(thruster_data(j).phi - ...
                dt * thruster_data(j).dphi_max, thruster_data(j).phi_min);
            numOfAzi = numOfAzi + 1;
    end
end
Nvar = numOfTunnel + 2 * numOfAzi + 3;%number of variables in quad
N_A = numOfTunnel * 2 + numOfAzi * 5;

H = eye(Nvar);
f= zeros(Nvar,1);
N_Aeq = 3;
Aeq= zeros(N_Aeq,Nvar);
beq=zeros(N_Aeq,1);
beq(1) = T_r.Tx;
beq(2) = T_r.Ty;
beq(3) = T_r.Tm;
A=zeros(N_A,Nvar);
b= zeros(N_A,1);

i_H = 1;
numOfTunnel=0;
numOfAzi=0;
pos_A = 1;
for j=1:N_enabled_thruster
    switch thruster_data(j).type
        case 4  %tunnel
            H(i_H,i_H) = thruster_data(j).weight;
            i_H = i_H + 1;
            
            for i = 1:3 %R1 R2 R3
                switch i
                    case 1
                        Aeq(i,numOfAzi * 2 + numOfTunnel + 1) = 0;
                    case 2
                        Aeq(i,numOfAzi * 2 + numOfTunnel + 1) = 1;
                    case 3
                        Aeq(i,numOfAzi * 2 + numOfTunnel + 1) = ...
                            thruster_data(j).x - thruster_data(j).x0;         
                end
            end
            
            A(pos_A,numOfAzi * 2 + numOfTunnel + 1) = -1;
            b(pos_A) = -thruster_data(j).Tplus;
            A(pos_A,numOfAzi * 2 + numOfTunnel + 1) = 1;
            b(pos_A) = thruster_data(j).T_;
            
            pos_A = pos_A + 2;
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
                                Aeq(i,numOfAzi * 2 + numOfTunnel + ixy) = 1;
                            case 2
                                Aeq(i,numOfAzi * 2 + numOfTunnel + ixy) = 0;
                            case 3
                                Aeq(i,numOfAzi * 2 + numOfTunnel + ixy) = ...
                                    - thruster_data(j).y + thruster_data(j).y0;

                        end
                        
                    elseif ixy ==2
                        switch i
                            case 1
                                Aeq(i,numOfAzi * 2 + numOfTunnel + ixy) = 0;
                            case 2
                                Aeq(i,numOfAzi * 2 + numOfTunnel + ixy) = 1;
                            case 3
                                Aeq(i,numOfAzi * 2 + numOfTunnel + ixy) = ...
                                    thruster_data(j).x - thruster_data(j).x0;

                        end                        
                    end
                end
                
                if ixy == 1
                    A(pos_A,numOfAzi * 2 + numOfTunnel + ixy) = -sin(thruster_data(j).phi_);%(3.15a)
                    b(pos_A) = -c2 * abs(cos(thruster_data(j).phi_)) * ...
                        (thruster_data(j).Tplus - thruster_data(j).T_);
                    
                    A(pos_A + 1,numOfAzi * 2 + numOfTunnel + ixy) = sin(thruster_data(j).phiPlus);%(3.15b)
                    b(pos_A + 1) = -c2 * abs(cos(thruster_data(j).phiPlus)) * ...
                        (thruster_data(j).Tplus - thruster_data(j).T_);
                    
                    A(pos_A + 2,numOfAzi * 2 + numOfTunnel + ixy) = cos(thruster_data(j).phi);%(3.16a)
                    b(pos_A + 2) = thruster_data(j).T_;
                    
                    A(pos_A + 3,numOfAzi * 2 + numOfTunnel + ixy) = -cos(thruster_data(j).phi_);%(3.16b)
                    b(pos_A + 3) = -cos(thruster_data(j).phi_ - thruster_data(j).phi) * ...
                        thruster_data(j).Tplus;
                     
                    A(pos_A + 4,numOfAzi * 2 + numOfTunnel + ixy) = -cos(thruster_data(j).phiPlus);%(3.16c)
                    b(pos_A + 4) = -cos(thruster_data(j).phiPlus - thruster_data(j).phi) * ...
                        thruster_data(j).Tplus;
                    
                elseif ixy == 2
                    A(pos_A,numOfAzi * 2 + numOfTunnel + ixy) = cos(thruster_data(j).phi_);%(3.15a)
                                  
                    A(pos_A + 1,numOfAzi * 2 + numOfTunnel + ixy) = -cos(thruster_data(j).phiPlus);%(3.15b)
                   
                    A(pos_A + 2,numOfAzi * 2 + numOfTunnel + ixy) = sin(thruster_data(j).phi);%(3.16a)
                    
                    A(pos_A + 3,numOfAzi * 2 + numOfTunnel + ixy) = -sin(thruster_data(j).phi_);%(3.16b)
                    
                    A(pos_A + 4,numOfAzi * 2 + numOfTunnel + ixy) = -sin(thruster_data(j).phiPlus);%(3.16c)
                   
                end
                
            end
            pos_A = pos_A + 5;

            

            numOfAzi = numOfAzi + 1;
    end
end

Aeq(1,numOfTunnel + 2 * numOfAzi + 1) = -1;
Aeq(2,numOfTunnel + 2 * numOfAzi + 2) = -1;
Aeq(3,numOfTunnel + 2 * numOfAzi + 3) = -1;

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
[solution,status] = Mpcqpsolver(Linv,f,A,b,Aeq,beq,iA0,opt);
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
            alloc_out(i).Ty = solution(numOfTunnel + numOfAzi * 2 + 1);
            alloc_out(i).T = alloc_out(i).Ty;
            alloc_out(i).Tm = alloc_out(i).Ty * (thruster_data(i).x ...
                - thruster_data(i). x0) + alloc_out(i).Tx ...
                * (thruster_data(i).y0 - thruster_data(i).y);
            
            numOfTunnel = numOfTunnel + 1;
        case 3
            alloc_out(i).Tx = solution(numOfTunnel + numOfAzi * 2 + 1);
            alloc_out(i).Ty = solution(numOfTunnel + numOfAzi * 2 + 2);
            alloc_out(i).T = sqrt(alloc_out(i).Tx^2 + alloc_out(i).Ty^2);
            alloc_out(i).phi = atan2(solution(numOfTunnel + numOfAzi * 2 + 2), ...
                solution(numOfTunnel + numOfAzi * 2 + 1));
            alloc_out(i).Tm = alloc_out(i).Ty * (thruster_data(i).x ...
                - thruster_data(i). x0) + alloc_out(i).Tx ...
                * (thruster_data(i).y0 - thruster_data(i).y);
           numOfAzi = numOfAzi + 1; 
    end
            
            
end

  