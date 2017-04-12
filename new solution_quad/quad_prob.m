[x0_1, x0_2, x0_3, x0_4 ,x0_5, x0_6 ]=deal(1.243,1.916,1.822,18.29,8.35,1.69);
x0=[x0_1, x0_2, x0_3, x0_4 ,x0_5, x0_6 ]';
[w1, w2, w3, w4 ,w5, w6 ]=deal(100,50,10,5,1,0.5);
w=[w1, w2, w3, w4 ,w5, w6 ]';
Q=diag(w);

[n_1,n_2]=deal(14.685850521606639,15.832079691957729);
n0=[n_1,n_2]';
n1=n0+[0.05 -0.02]';
D_n=n1-n0;
h=[constr1(x0,n_1)
    constr2(x0,n_2)];

% Jh=subs(Jh,x,x0);
% A=[Q Jh'
%     Jh zeros(size(Jh,1),size(Jh,1))];
% d=[0.01;0.042];
% B=[-c;d];
% X=vpa(A\B,5)
c=zeros(size(x0));
iter=100;
d=D_n/iter;

V_x0=zeros(6,1);
dt=0.2;
[UB1, UB2, UB3, UB4, UB5, UB6]=deal(pi,pi,pi,21.89,10.709,1.907);
UB=[UB1, UB2, UB3, UB4, UB5, UB6]';
[LB1 LB2 LB3 LB4 LB5 LB6]=deal(pi/4,pi/4,pi/2,21.89-6.39,6.709,1.197);
LB=[LB1 LB2 LB3 LB4 LB5 LB6]';
[V_limt1, V_limt2, V_limt3, V_limt4, V_limt5, V_limt6]=deal(5/180*pi, 5/180*pi, 5/180*pi, 0.1, 0.1, 0.1)
V_limt=[V_limt1, V_limt2, V_limt3, V_limt4, V_limt5, V_limt6]';
[A_limt1, A_limt2, A_limt3, A_limt4, A_limt5, A_limt6]=deal(5/180*pi, 5/180*pi, 5/180*pi, 0.1, 0.1, 0.1);
A_limt=[A_limt1, A_limt2, A_limt3, A_limt4, A_limt5, A_limt6]';
UB_v=V_limt*dt+x0;
LB_v=-V_limt*dt+x0;
UB_a=A_limt*dt*dt+V_x0*dt+x0;
LB_a=-A_limt*dt*dt+V_x0*dt+x0;


for i=1:iter
    i
    Q=diag(w);
    J=Jh(x0);
    c=zeros(size(x0));
    d_x=zeros(size(x0));
    %     d_xlim=zeros(size(x0));%for speed limit and acc limit
    %     d=d-d_xlim*x;
    %     cut_x=[2 4];
    %     keep_x=[1 3 5 6];

    
    A=double([Q J'
        J zeros(size(J,1),size(J,1))]);
    
    B=[-c;d];
    temp=A\B;
    temp(end-1:end)=[];
    d_x=temp;
    
    x0_=x0;
    x0=x0_+d_x;
 %% 
LB_v=-V_limt*dt+x0;
UB_a=A_limt*dt*dt+V_x0*dt+x0;
LB_a=-A_limt*dt*dt+V_x0*dt+x0;

 cut_x=zeros(size(x0));
 keep_x=ones(size(x0));

while 1
    temp=(x0<min([LB LB_v LB_a],[],2) | x0>max([UB UB_v UB_a],[],2));
    if max(temp)==0
        break;%solution is in valid range
    end
    cut_x=cut_x|temp;
    keep_x=~cut_x;
    cut_num=cut_x.*(1:6)';
    keep_num=keep_x.*(1:6)';
    flag_cutx=(max(cut_x)~=0);
    Q(cut_num(cut_num>0),:)=[];
    Q(:,cut_num(cut_num>0))=[];
    c(cut_num(cut_num>0),:)=[];
    J(:,cut_num(cut_num>0))=[];
        A=double([Q J'
        J zeros(size(J,1),size(J,1))]);
    
    B=[-c;d];
    temp=A\B;
    temp(end-1:end)=[];
    d_x(keep_num(keep_num>0),1)=temp;
    x0=x0_+d_x;
end
 %%   
    
    h=[constr1(x0,0)
        constr2(x0,0)];
    if i~=iter
        d=(n1-h)/(iter-i);
    end
end
h=[constr1(x0,n1(1))
    constr2(x0,n1(2))];