points=[0,0;53124.5000000000,-24810.5000000000;55055,-24667.5000000000;57986.5000000000,-24381.5000000000;60060,-21450;63992.5000000000,-20520.5000000000;65780,-14300;68354,-11654.5000000000;68997.5000000000,-5720;70070,0;68997.5000000000,5720;68354,11654.5000000000;65780,14300;63992.5000000000,20520.5000000000;60060,21450;57986.5000000000,24381.5000000000;55055,24667.5000000000;53124.5000000000,24810.5000000000;0,0];

[max_dat]=max(points,[],1);
[min_dat]=min(points,[],1);

[x y]=meshgrid(min_dat(1):1000:max_dat(1),min_dat(2):1000:max_dat(2));
M=size(x,1);
NN=size(x,2);
dat_sav = [];

N = size(points,1) - 1;
b = zeros(N,1);
A = zeros(N,2);
%     points(end + 1,:) = points(1,:);
for i = 1 : N
    A(i,1) = -points(i + 1,2) + points(i,2);
    A(i,2) = -points(i,1) + points(i + 1,1);
    b(i) = -points(i,1) * points(i + 1,2) + points(i + 1,1) * points(i,2);


end

plot(points(:,2),points(:,1));
axis equal
hold on
for i = 1:M
    for j=1:NN
    if A * [x(i,j) y(i,j)]' >= b
        dat_sav(end+1,:) = [x(i,j) y(i,j)];
        plot(y(i,j),x(i,j),'o');
    end
    end
end
hold off