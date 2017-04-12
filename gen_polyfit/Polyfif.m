function [r,f]=Polyfif(x,y,N_poly,x_new)
coder.inline('never')
 [r, ~] = polyfit(x, y,N_poly);
 
 f = polyval(r,x_new);
 coder.inline('never')
end