function y=inf_test(x)
persistent i;
if isempty(i)
    i = 0;
end


i = i + 1;
if i >= 100
    y = 1;
else
    y = sqrt(x + inf_test(x));
end
end