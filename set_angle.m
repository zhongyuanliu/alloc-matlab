function y = set_angle(thruster_data_in,alloc_out)
%% set angle for azi
temp = thruster_data_in.dphi_max * thruster_data_in.dt;
if angle_dist(thruster_data_in.phi,alloc_out.phi) >= ...
        angle_dist(thruster_data_in.phi,thruster_data_in.phi_min) ...
        && angle_dist(thruster_data_in.phi,alloc_out.phi) ...
        >=angle_dist(thruster_data_in.phi,thruster_data_in.phi_max)
    
    if 2 * pi - angle_dist(thruster_data_in.phi,alloc_out.phi) < temp
        y = alloc_out.phi;
    else
        y = Pi_toPi(thruster_data_in.phi - temp);
    end
else
    if angle_dist(thruster_data_in.phi,alloc_out.phi) < temp
        y = alloc_out.phi;
    else
        y = Pi_toPi(thruster_data_in.phi + temp);
    end    
    
end

end