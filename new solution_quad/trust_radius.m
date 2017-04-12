function y=trust_radius(Delta,Delta_min,Delta_max,s,ared,pred)
    r=ared/pred;
    if r<1e-4 || ared<0 
        beta1=0.5*(Delta<=10)+0.25*(Delta<=10e4&&Delta>10)+0.1*(Delta>1e4);
        Delta=max([beta1*min(Delta,norm(s)) Delta_min]);
    else
        
    end
    if r>=1e-4 && r<1e-2
        Delta=max([0.5*Delta Delta_min norm(s)]);
    end
    if r>=1e-2 && r<0.1
        Delta=max([0.9*Delta Delta_min]);
    end
    if r>=0.75
        beta2=10*(Delta<=1e-4)+5*(Delta>1e-4&&Delta<=0.1)...
            +4*(Delta>0.1&&Delta<=100)+2*(Delta>100);
        Delta=min([beta2*Delta Delta_max]);
    end    
    y=Delta;
end