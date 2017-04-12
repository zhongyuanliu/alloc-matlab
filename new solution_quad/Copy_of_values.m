%%
clear
[w1, w2, w3, w4, w5, w6 ]=deal(1000,1000,100,100,10,10);
[x1, x2, x3, x4, x5, x6 ]=deal(1.243,1.916,1.822,18.29,8.35,1.69);
x=[x1, x2, x3, x4, x5, x6 ]';
[x0_1, x0_2, x0_3, x0_4 ,x0_5, x0_6 ]=deal(1.243,1.916,1.822,18.29,8.35,1.69);
x0=[x0_1, x0_2, x0_3, x0_4 ,x0_5, x0_6 ]';
v=1'
[ V_x0_1 ,V_x0_2 ,V_x0_3, V_x0_4 ,V_x0_5, V_x0_6 ]=deal(0);
[LB1, LB2, LB3, LB4, LB5, LB6 ]=deal(pi/4,pi/4,pi/2,21.89-6.39,6.709,1.19);
[UB1, UB2, UB3, UB4, UB5, UB6]=deal(pi,pi,pi,21.89,10.709,1.907);
UB=[UB1, UB2, UB3, UB4, UB5, UB6]';
[LB1 LB2 LB3 LB4 LB5 LB6]=deal(pi/4,pi/4,pi/2,21.89-6.39,6.709,1.197);
LB=[LB1 LB2 LB3 LB4 LB5 LB6]';
[V_limt1, V_limt2, V_limt3, V_limt4, V_limt5, V_limt6]=deal(5/180*pi, 5/180*pi, 5/180*pi, 0.1, 0.1, 0.1)
V_limt=[V_limt1, V_limt2, V_limt3, V_limt4, V_limt5, V_limt6]';
[A_limt1, A_limt2, A_limt3, A_limt4, A_limt5, A_limt6]=deal(5/180*pi, 5/180*pi, 5/180*pi, 0.1, 0.1, 0.1);
A_limt=[A_limt1, A_limt2, A_limt3, A_limt4, A_limt5, A_limt6]';
 [n_1,n_2]=deal(14.687,15.85);
 n=[n_1,n_2]';
 V_x0=zeros(6,1);
 A_x0=zeros(6,1);
 miu=1;%???????????????
 
[s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 s17...
    s18 s19 s20 s21 s22 s23 s24 s25 s26 s27 s28 s29 s30 s31...
    s32 s33 s34 s35 s36]=deal(0.02);
s=[s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 s17...
    s18 s19 s20 s21 s22 s23 s24 s25 s26 s27 s28 s29 s30 s31...
    s32 s33 s34 s35 s36]';
 [y1, y2 ]=deal(1);
 lambda_h= [y1, y2 ]';
[ lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14 lambda15 lambda16 lambda17...
    lambda18 lambda19 lambda20 lambda21 lambda22 lambda23 lambda24 lambda25 lambda26 lambda27 lambda28 lambda29 lambda30 lambda31...
    lambda32 lambda33 lambda34 lambda35 lambda36]=deal(0.6);
lambda_g=[ lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14 lambda15 lambda16 lambda17...
    lambda18 lambda19 lambda20 lambda21 lambda22 lambda23 lambda24 lambda25 lambda26 lambda27 lambda28 lambda29 lambda30 lambda31...
    lambda32 lambda33 lambda34 lambda35 lambda36]';
%%
