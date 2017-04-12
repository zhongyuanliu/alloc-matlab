function y=test_gen(x)
% A=zeros(3,1);
% if x>0
%     A(end)=x
% else
%     if ~isempty(A)
%     A(1)=[];
%     end
% end
% y=A;
coder.extrinsic('mat2cell')
[m,n] = size(x);
y=cell(1);
y{1}=1;
y{end+1}=1;
y=mat2cell(x,1:m,1:n)

end