Copy_of_values;
size_z=42;
size_lambda=38;
size_x=6;
size_s=36;
size_lambda_h=2;
size_lambda_g=36;
v=1;
z=[x;s];
lambda=[lambda_h;lambda_g];

tau=0.995;%tau is between 0 and 1, typically close to one

% [y,Dy]=merit_fcn(x,s,d_z,miu,x0,n,v,UB,LB,V_limt,A_limt,V_x0)
Delta=0.01;
delta=1e-5;%threshold
Delta_max=100;%not sure
Delta_min=delta;%not sure
eta=1e-8;
imax=3;
alpha_T=0.5;




miu=1;
rho=0.1;%
[AA,BB,A,B,LB,UB,A_x,W]=values(x,s,lambda_h,lambda_g,miu);
h=[constr1(x,n_1)
    constr2(x,n_2)];
g=[x-UB
    LB-x
    A*x-B];
k=0;
while (sum(abs(h)>1e-7)+sum(g>0))>=1
    h=[constr1(x,n_1)
        constr2(x,n_2)];
    g=[x-UB
        LB-x
        A*x-B];
    
    while (sum(abs(h)>1e-7)+sum(abs(g-s)>1e-7))>=1
        h=[constr1(x,n_1)
            constr2(x,n_2)];
        g=[x-UB
            LB-x
            A*x-B];
        
        
        x=z(1:size_x,1);
        s=z(size_x+1:size_x+size_s,1);
        [AA,BB,A,B,LB,UB,A_x,W]=values(x,s,lambda_h,lambda_g,miu);
        neig=sum(eig(AA)<0);
        LineSearch=0;
        
        if neig<=size(h,1)+size(g,1)
            d=AA\BB;
            d_s=d(size_x+1:size_x+size_s,1);
            d_z=d(1:size_z,1);
            d_x=d(1:size_x,1);
            d_lambda=d(size_z+1:end,1);
            d_h=d(size_z+1:size_z+size_lambda_h,1);
            d_g=d(size_z+size_lambda_h+1:end,1);
            alpha_z_max=min((d_s<0).*(-tau*s./d_s)+(d_s>=0));
            alpha_lambda_max=min((d_g<0).*(-tau*lambda_g./d_g)+(d_g>=0));
            
%             if min([alpha_z_max alpha_lambda_max])>delta
                j=0;
                alpha_T=1;
                
                while(j<=imax && alpha_T>delta && LineSearch==0 )
                    [y,Dy,c]=merit_fcn(z,d_z,miu,x0,n,v,UB,LB,V_limt,A_limt,V_x0);
                    [y1,Dy1,c1]=merit_fcn(z+alpha_T*alpha_z_max*d_z,d_z,miu,x0,n,v,UB,LB,V_limt,A_limt,V_x0);
                    alpha_T=min(0.5,Delta/norm(d_z));
                    alpha=(d_z'*W*d_z>0);
                    v_trial=(grad_fmiu(x,s,x0,miu)'*d_z+alpha/2*d_z'*W*d_z)/((1-rho)*norm(c));
                    v=v*(v>=v_trial)+(v_trial+1)*(v<v_trial);
                    
                    
                    if (y1<=(y+eta*alpha_T*alpha_z_max*Dy))
                        alpha_z=alpha_T*alpha_z_max;
                        alpha_lambda=alpha_T*alpha_lambda_max;
                        
                        
%                         ared_dz=merit_fcn(z,d_z,miu,x0,n,v,UB,LB,V_limt,A_limt,V_x0)...
%                             -merit_fcn(z+d_z,d_z,miu,x0,n,v,UB,LB,V_limt,A_limt,V_x0);
                        ared_dz=merit_fcn(z,alpha_z*d_z,miu,x0,n,v,UB,LB,V_limt,A_limt,V_x0)...
                            -merit_fcn(z+alpha_z*d_z,alpha_z*d_z,miu,x0,n,v,UB,LB,V_limt,A_limt,V_x0);
          
                        pred_dz=-grad_fmiu(x,s,x0,miu)'*alpha_z*d_z-alpha/2*alpha_z*d_z'*W*alpha_z*d_z+v*(norm(c)-norm(c+A_x*alpha_z*d_z));
                       
                        z=z+alpha_z*d_z;
                        x=z(1:size_x,1)
                        s=z(size_x+1:size_x+size_s,1);
                        lambda=lambda+alpha_lambda*d_lambda;
                        lambda_h=lambda(1:size_lambda_h,1);
                        lambda_g=lambda(size_lambda_h+1:end,1);
                        Delta=trust_radius(Delta,Delta_min,Delta_max,s,ared_dz,pred_dz)
                        LineSearch=1;
                    else
                        j=j+1;
                        alpha_T=min(0.5,Delta/norm(d_z));
                    end
                    
                end
%             end
        end
        if LineSearch==0
            %             Delta=2*alpha_z*norm(d_z);
        end
        k=k+1;
    end
    miu=miu/100*(j<=3)+miu/5*(j>3);%??
    
end