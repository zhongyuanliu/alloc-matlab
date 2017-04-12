clear
def_struct
clear qpsolver
% tic
% T = [426.277880858579;286.678460356666;446.239669612964;266.716671601555;0;346.497276230906;275.034083582627;0]';
% phi = [-0.0478344712504149;0.0500000000000000;0.0500000000000000;-0.0499999999999998;0;-0.0499999999999998;0.0500000000000000;0]';
TT = 0;
set_T = struct('T',{0,0,0,0,0,0,0,0},'phi',0);
T_r.Tx = 10000;
T_r.Ty = 0;
%     T_r = struct('Tx',100000,'Ty',0,'Tm',100000);
T_r.Tm = 0;
k=12;
T_r.Tx = T_r.Tx * k;
T_r.Ty = T_r.Ty * k;
T_r.Tm = T_r.Tm * k;
T4 = [];
init = true;
    thruster_data_in = thruster_data0;
for j = 1:1
    N_enabled_thruster=0;
        T_r.Tx = 3e4 ; %sin((j-1)/20)*1;
        T_r.Ty = -13e4*1 ;
%         T_r.Tx = 2.8e3;
%         T_r.Ty = -12.7e3;
        T_r.Tm = -146.7e3*0;
    % %     T_r = struct('Tx',100000,'Ty',0,'Tm',100000);
    
    
    

    for i = 1:8
        %     thruster_data(i).enable = round(rand(1));
        if thruster_data(i).enable == 1
            N_enabled_thruster = N_enabled_thruster + 1;
            if j == 1
             thruster_data_in(N_enabled_thruster) = thruster_data(i);
            end
            %         thruster_data_in(N_enabled_thruster).T = T(i);
            %         thruster_data_in(N_enabled_thruster).phi = phi(i);
            if j>1
                %             set_T(i).T = alloc_out(i).T
%                 set_T(N_enabled_thruster).phi = set_angle(thruster_data_in(N_enabled_thruster), alloc_out(N_enabled_thruster));
                thruster_data_in(N_enabled_thruster).T = alloc_out(N_enabled_thruster).T;
                thruster_data_in(N_enabled_thruster).phi =  alloc_out(N_enabled_thruster).phi;
                thruster_data_in(N_enabled_thruster).T_reserve = alloc_out(N_enabled_thruster).T;
                thruster_data_in(N_enabled_thruster).phi_reserve =  alloc_out(N_enabled_thruster).phi;
                %             thruster_data_in(N_enabled_thruster).phi =  set_T(N_enabled_thruster).phi;
            end
        end
    end
    
    [solution,alloc_out,thruster_data_in,status]=qpsolver(init,thruster_data_in,T_r,N_enabled_thruster,rudder_table,0,1);
    init = false;
    if length(solution)>2
        error = (norm(solution(end-2:end)))
    end
    [alloc_out.T];
    % [alloc_out.phi]
    status
    T4(end+1) = alloc_out(4).T;
end
sum([alloc_out.Tm])
% toc
% figure;
% plot(1:numel(T4),T4)
