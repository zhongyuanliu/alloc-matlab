function [solution,status,iA,lambda]= Mpcqpsolver(Linv,f,A,b,Aeq,beq,iA0,opt)
coder.inline('never');
[solution,status,iA,lambda] = mpcqpsolver(Linv,f,A,b,Aeq,beq,iA0,opt);
end
