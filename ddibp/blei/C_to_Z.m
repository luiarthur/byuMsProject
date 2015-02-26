function z = C_to_Z(c,kn)
    
    z = zeros(size(c));
    for n = 1:length(c)
        if reachable(c,n,kn)
            z(n) = 1;
        else
            z(n) = 0;
        end
    end
    
    z(kn)=1;