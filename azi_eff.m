function y = azi_eff(Ce,phi)
% y = 0;
% N = size(Ce,1);
phi = to2Pi(phi);

% if phi>=Ce(end-1,1)
%     y1 = Ce(end-1,2);
%     y2 = Ce(end,2);
% 
%     y = (phi - Ce(end-1,1)) / (Ce(end,1) - Ce(end-1,1)) * (y2 - y1) + y1;
% else
%     t= ((phi >= Ce(1:end-1,1) & phi<Ce(2:end,1)));
%     I = find(t==1);
%     y1 = Ce(I,2);
%     y2 = Ce(I + 1,2);
%     y = (phi - Ce(I,1)) / (Ce(I + 1,1) - Ce(I,1)) * (y2 - y1) + y1;
% end
% % % % % % % % % % % % % % % % % % % % % % % % % 
if phi>=Ce(end-1,1)
    y = Ce(end-1,2);
else
    t= ((phi >= Ce(1:end-1,1) & phi<Ce(2:end,1)));
    y = Ce(t==1,2);
end

% for i = 1 : N - 1
%     if (phi>= Ce(i,1) && phi<Ce(i + 1,1)) ...
%         || phi>=Ce(end-1,1)
%         y = Ce(i,2);
%     end
% end


end