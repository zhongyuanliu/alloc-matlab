clear
syms x1 x2 x3 x4 x5 x6 real
syms x0_1 x0_2 x0_3 x0_4 x0_5 x0_6 real
syms V_x0_1 V_x0_2 V_x0_3 V_x0_4 V_x0_5 V_x0_6 real
syms LB1 LB2 LB3 LB4 LB5 LB6 real
syms UB1 UB2 UB3 UB4 UB5 UB6 real
syms V_limt1 V_limt2 V_limt3 V_limt4 V_limt5 V_limt6 real
syms A_limt1 A_limt2 A_limt3 A_limt4 A_limt5 A_limt6 real
syms n_1 n_2 miu real
syms s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 s17...
    s18 s19 s20 s21 s22 s23 s24 s25 s26 s27 s28 s29 s30 s31...
    s32 s33 s34 s35 s36 real
syms y1 y2 real
syms lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14 lambda15 lambda16 lambda17...
    lambda18 lambda19 lambda20 lambda21 lambda22 lambda23 lambda24 lambda25 lambda26 lambda27 lambda28 lambda29 lambda30 lambda31...
    lambda32 lambda33 lambda34 lambda35 lambda36  real
syms n miu real
syms w1 w2 w3 w4 w5 w6 real
%lambda denotes the Lagrange multiplier vector associated with constraints g
%y denotes the Lagrange multiplier vector associated with h.

s=[s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 s17...
    s18 s19 s20 s21 s22 s23 s24 s25 s26 s27 s28 s29 s30 s31...
    s32 s33 s34 s35 s36 ]'
y=[y1 y2]';
lambda=[lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14 lambda15 lambda16 lambda17...
    lambda18 lambda19 lambda20 lambda21 lambda22 lambda23 lambda24 lambda25 lambda26 lambda27 lambda28 lambda29 lambda30 lambda31...
    lambda32 lambda33 lambda34 lambda35 lambda36 ]';
x=[x1,x2,x3,x4,x5,x6]'
w=[w1 w2 w3 w4 w5 w6]';
x0=[x0_1 x0_2 x0_3 x0_4 x0_5 x0_6]'

V_x0=[V_x0_1 V_x0_2 V_x0_3 V_x0_4 V_x0_5 V_x0_6]'
LB =[LB1 LB2 LB3 LB4 LB5 LB6]'
UB =[UB1 UB2 UB3 UB4 UB5 UB6]'
V_limt =[V_limt1 V_limt2 V_limt3 V_limt4 V_limt5 V_limt6]'
A_limt =[A_limt1 A_limt2 A_limt3 A_limt4 A_limt5 A_limt6]'


%% linear contr
dt=0.2
% x0=[pi;pi;pi;21.89;10.709;1.907];
% V_x0=zeros(6,1);
% LB=[pi/4;pi/4;pi/2;21.89-6.39;6.709;1.197];
% UB=[pi;pi;pi;21.89;10.709;1.907];
% V_limt=[5/180*pi 5/180*pi 5/180*pi 0.1 0.1 0.1]';
% A_limt=[5/180*pi 5/180*pi 5/180*pi 0.1 0.1 0.1]'*2;
A=[eye(6);-eye(6);eye(6);-eye(6)];
% -[5/180*pi 5/180*pi 5/180*pi 0.1 0.1 0.1]'*dt<=x-x0<=[5/180*pi 5/180*pi 5/180*pi 0.1 0.1 0.1]'*dt;
B=[V_limt*dt+x0
    V_limt*dt-x0 %speed restraint
    A_limt*dt*dt+V_x0*dt+x0 %accel restraint
    A_limt*dt*dt-V_x0*dt-x0
    ];

UB_v=V_limt*dt+x0;
LB_v=-V_limt*dt+x0;
UB_a=A_limt*dt*dt+V_x0*dt+x0;
LB_a=-A_limt*dt*dt+V_x0*dt+x0;
x_L=[LB;LB_v;LB_a];
x_U=[UB;UB_v;UB_a];
%%
% costfun(x1,x2,x3,x4,x5,x6,x0_1,x0_2,x0_3,x0_4,x0_5,x0_6)+constr1(x1,x2,x3,x4,x5,x6,n)
% A*[x1 x2 x3 x4 x5 x6 ]'-B
h=[constr1(x,n_1)
    constr2(x,n_2)];

g=[x-UB
    LB-x
    A*x-B];
f_cost=costfun(x,x0,w);
grad_cost=gradient(f_cost,x);
% s=-g;
% gradient(h(1),x);
% gradient(h(2),x);
% hessian(h(1),x);
% hessian(h(2),x);
syms obj_factor real
Jh=jacobian(h,x);
Jg=jacobian(g,x);

f=costfun(x,x0,w);
f_lagr=obj_factor*f+lambda1*h(1)+lambda2*h(2);
H_lagr=hessian(f_lagr,x)
% fmiu=costfun(x,x0)-miu*sum(log(s));



















% f_hss=hessian(f,x);
% temp=0;
% 
% % dd_g=sym('dd_g',[size(x,1) size(g,1)]);
% g_hss=sym(zeros(size(x,1),size(x,1), size(g,1)));
% G_hss=sym(zeros(size(x,1),size(x,1)));
% for i=1:size(g,1)
%     %     for j=1:size(x,1)
%     g_hss(:,:,i)=lambda(i)*hessian(g(i),x);
%     G_hss=g_hss(:,:,i)+G_hss;
%     %     end
% end
% 
% 
% h_ss=sym(zeros(size(x,1),size(x,1), size(h,1)));
% H_ss=sym(zeros(size(x,1),size(x,1)));
% for i=1:size(h,1)
%     %     for j=1:size(x,1)
%     h_ss(:,:,i)=y(i)*hessian(h(i),x);
%     H_ss=H_ss+h_ss(:,:,i);
%     %     end
% end
% 
% H=f_hss+G_hss+H_ss;%check 
% % for i=1:size(g,1)
% % temp=lambda(i)*diff(g,x(i),2)+temp;
% % end
% %%
% S=diag(s);
% Lambda=diag(lambda);
% % AA=[H sym(zeros(size(H,1),size(S,2))) Jh' Jg'
% %     sym(zeros(size(S,1),size(H,2))) S*Lambda sym(zeros(size(S,1),size(Jh,1))) -S
% %     Jh sym(zeros(size(Jh,1),size(S,2))) sym(eye(size(Jh,1),size(Jh,1))) sym(zeros(size(Jh,1),size(S,2)))
% %     Jg -S sym(zeros(size(Jg,1),size(Jh,1))) sym(eye(size(S,1),size(S,2)))];
% W=[H sym(zeros(size(H,1),size(S,2)))
%    sym(zeros(size(S,1),size(H,2))) (S)\Lambda];
% A_x=[Jh sym(zeros(size(Jh,1),size(S,2)))
%     Jg sym(eye(size(Jg,1),size(S,2)))];
% AA=[H sym(zeros(size(H,1),size(S,2))) Jh' Jg'
%     sym(zeros(size(S,1),size(H,2))) Lambda sym(zeros(size(S,1),size(Jh,1))) S
%     Jh sym(zeros(size(Jh,1),size(S,2))) sym(zeros(size(Jh,1),size(Jh,1))) sym(zeros(size(Jh,1),size(Jg,1)))
%     Jg sym(eye(size(Jg,1),size(S,2))) sym(zeros(size(Jg,1),size(Jh,1))) sym(zeros(size(S,1),size(Jg,1)))];
% 
% BB=-[gradient(f,x)+Jh'*y+Jg'*lambda
%     S*lambda-miu*ones(size(g))
%     h
%     g+s]
% size_x=size(H,1)
% size_s=size(S,1)
% size_z=size(H,1)+size(S,1);
% size_lambda=size(Jh,1)+size(Jg,1);
% size_h=size(Jh,1);
% size_g=size(Jg,1);

% %% Values
% [x1, x2, x3, x4, x5, x6 ]=deal(pi/3,2,1.9,25,5,1);
% [x0_1, x0_2, x0_3, x0_4 ,x0_5, x0_6 ]=deal(pi/3,2,1.9,25,5,1);
% [ V_x0_1 ,V_x0_2 ,V_x0_3, V_x0_4 ,V_x0_5, V_x0_6 ]=deal(0);
% [LB1, LB2, LB3, LB4, LB5, LB6 ]=deal(pi/4,pi/4,pi/2,21.89-6.39,6.709,1.19);
% [UB1, UB2, UB3, UB4, UB5, UB6]=deal(pi,pi,pi,21.89,10.709,1.907);
% [V_limt1, V_limt2, V_limt3, V_limt4, V_limt5, V_limt6]=deal(5/180*pi, 5/180*pi, 5/180*pi, 0.1, 0.1, 0.1)
% [A_limt1, A_limt2, A_limt3, A_limt4, A_limt5, A_limt6]=deal(5/180*pi, 5/180*pi, 5/180*pi, 0.1, 0.1, 0.1);
%  n_1=[12];
%  n_2=18;
%  miu=1;%???????????????
%  
% [s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12 s13 s14 s15 s16 s17...
%     s18 s19 s20 s21 s22 s23 s24 s25 s26 s27 s28 s29 s30 s31...
%     s32 s33 s34 s35 s36]=deal(1.2);
%  [y1, y2 ]=deal(1);
% [ lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14 lambda15 lambda16 lambda17...
%     lambda18 lambda19 lambda20 lambda21 lambda22 lambda23 lambda24 lambda25 lambda26 lambda27 lambda28 lambda29 lambda30 lambda31...
%     lambda32 lambda33 lambda34 lambda35 lambda36]=deal(0.6)
% 
