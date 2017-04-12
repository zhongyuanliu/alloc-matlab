
%!!!!!!!!!!!!!check if there are bugs
function [phi_,phiPlus]=angleMaxMin(angle_start,angle_end,phi,phi1,phi2)
%% up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise
% phi: current angle
% phi1 to phi2
%!!!!!!!!!!!!!!!!!!!!!watch out for FeasibilityTol
FeasibilityTol = 1e-8;
angle_start = Pi_toPi(angle_start);
angle_end = Pi_toPi(angle_end);
phi = Pi_toPi(phi);
phi1 = Pi_toPi(phi1);
phi2 = Pi_toPi(phi2);
effective_dist = angle_dist(angle_start,angle_end);
start_phi1 = angle_dist(angle_start,phi1);
phi1_end = angle_dist(phi1,angle_end);
phi1_phi2 = angle_dist(phi1,phi2);
phi1_start = angle_dist(phi1,angle_start);
phi1_phi = angle_dist(phi1,phi);
start_phi = angle_dist(angle_start,phi);
start_phi2 = angle_dist(angle_start,phi2);
phi_ = 0;
phiPlus = 0;
if abs(angle_start - angle_end) < FeasibilityTol%angle_start == angle_end means no forbidden zone
    phi_ = phi1;
    phiPlus = phi2;
else
    
if effective_dist >= start_phi1 %phi1 is at the range from start to end
    if phi1_end >= phi1_phi2 %phi2 is at the range from phi1 to end
        phi_ = phi1;
        phiPlus = phi2;
    elseif phi1_phi2 > phi1_end ...
            && phi1_phi2 <= phi1_start
        phi_ = phi1;
        phiPlus = angle_end;
    else %phi2 is at the range from start to phi1
        if phi1_phi <= phi1_end%phi is at the range from phi1 to end
            phi_ = phi1;
            phiPlus = angle_end;
        elseif start_phi <= start_phi2%phi is at the range from start to phi2
            phi_ = angle_start;
            phiPlus = phi2;
        end
        
    end
else%phi1 is at the range from end to start
    phi_ = angle_start;
    phiPlus = phi2;
end

end


end
