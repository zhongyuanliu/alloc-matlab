function baseT_generate
T1 = [1160*2 1.04];
T2 = [1160*2 -1.04];
T3 = [1110*2 pi];
T1(2) = acos(T3(1) / 2 / T1(1));
T2(2) = -T1(2);
end