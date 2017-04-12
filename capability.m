function capability
def_struct
clear qpsolver

Fc = [-1.8e5 4.2e5 3.2e5]'*0;%current force
Fw0 = [1.8e3 0e3 0 %0
    1.6e3 1.2e3 3e3 %10
    1.3e3 2.5e3 6e3 %20
    1.0e3 4e3 1e4    %30
    0.95e3 4.5e3 3e4  %40
    0.92e3 5.2e3 10e4  %50
    0.8e3 6.3e3 13e4    %60
    0.7e3 7.2e3 12e4   %70
    0.4e3 7.9e3 10e4   %80
    0.1e3 8.8e3 3e4    %90
    -0.4e3 7.9e3 -10e4   %100
    -0.7e3 7.2e3 -12e4   %110
    -0.8e3 6.3e3 -13e4    %120
    -0.92e3 5.2e3 -10e4  %130
     -0.95e3 4.5e3 -3e4  %140
    -1.0e3 3.6e3 -1.2e4    %150
    -1.3e3 2.5e3 -6e3 %160
    -1.6e3 1.2e3 -3e3 %170
    -1.9e3 0e3 0    %180
    ];%basic wind force
Fw_ = Fw0;
Fw_(:,2) = -Fw_(:,2);
Fw_(:,3) = -Fw_(:,3);
Fw_(1,:) = [];
Fw_(end,:) = [];
Fw0 = [Fw0; fliplr(Fw_')'];
Fw0 = [(0:10:350)' Fw0];

% tic
% T = [426.277880858579;286.678460356666;446.239669612964;266.716671601555;0;346.497276230906;275.034083582627;0]';
% phi = [-0.0478344712504149;0.0500000000000000;0.0500000000000000;-0.0499999999999998;0;-0.0499999999999998;0.0500000000000000;0]';
TT = 70000;
set_T = struct('T',{0,0,0,0,0,0,0,0},'phi',0);
% T_r.Tx = -2000;
% T_r.Ty = 3000;
% %     T_r = struct('Tx',100000,'Ty',0,'Tm',100000);
% T_r.Tm = 0;

T4 = [];
init = true;
speed_cap = zeros(36,1);
for dir = 1:36
    
    k1=0;
    k2 = 50;
    kp0 = 0;
    kp = k2 / 2;
    j = 0;
    while 1
        
        j = j + 1;
        N_enabled_thruster=0;
        
        T_r.Tx = Fc(1) + kp^2 * Fw0(dir,2);
        T_r.Ty = Fc(2) + kp^2 * Fw0(dir,3);
        T_r.Tm = Fc(3) + kp^2 * Fw0(dir,4);
        
        thruster_data_in = thruster_data0;
        for i = 1:8
            %     thruster_data(i).enable = round(rand(1));
            if thruster_data(i).enable == 1
                N_enabled_thruster = N_enabled_thruster + 1;
                thruster_data_in(N_enabled_thruster) = thruster_data(i);
                %         thruster_data_in(N_enabled_thruster).T = T(i);
                %         thruster_data_in(N_enabled_thruster).phi = phi(i);
                if j>1
                    %             set_T(i).T = alloc_out(i).T
                    set_T(N_enabled_thruster).phi = set_angle(thruster_data_in(N_enabled_thruster), alloc_out(N_enabled_thruster));
                    thruster_data_in(N_enabled_thruster).T = alloc_out(N_enabled_thruster).T;
                    thruster_data_in(N_enabled_thruster).phi =  alloc_out(N_enabled_thruster).phi;
                    %             thruster_data_in(N_enabled_thruster).phi =  set_T(N_enabled_thruster).phi;
                end
            end
        end
        
        [solution,alloc_out,status]=qpsolver(init,thruster_data_in,T_r,N_enabled_thruster,rudder_table,0,1);
        init = false;
        
        
        
        if length(solution)>2
            error = sqrt(norm(solution(end-2:end)));
            if error < 500
                k1 = kp;
                
            else
                k2 = kp;
            end
            kp = (k1 + k2) / 2;
            
        end
        
        if abs(kp0 - kp)<0.5
            
            break;
        end
        kp0 = kp;%previous kp
        
        [alloc_out.T];
        % [alloc_out.phi]
        %     status
        %     T4(end+1) = alloc_out(4).T;
    end
    speed_cap(dir) = kp;
end

polarplot(Fw0(:,1) / 180 * pi,speed_cap)
end