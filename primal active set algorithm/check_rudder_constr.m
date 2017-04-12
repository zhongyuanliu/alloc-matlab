function check = check_rudder_constr(points,A,b)
%%
check = 0;
[max_dat]=max(points,[],1);
[min_dat]=min(points,[],1);
interv = (max_dat - min_dat)/40;
[x y]=meshgrid(min_dat(1):interv:max_dat(1),min_dat(2):interv:max_dat(2));
M=size(x,1);
NN=size(x,2);
dat_sav = [];

plot(points(:,2),points(:,1));
axis equal
hold on
for i = 1:M
    for j=1:NN
        if A * [x(i,j) y(i,j)]' < b
            dat_sav(end+1,:) = [x(i,j) y(i,j)];
            plot(x(i,j),y(i,j),'*');
            check = 1;
        end
    end
end
% hold off
end