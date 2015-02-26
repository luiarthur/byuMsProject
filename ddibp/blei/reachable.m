function r = reachable(c,i,j)
    %determines whether j is reachable from i in c (column of C)
    
    q = i; v = i; stop = 0; r = 0;
    
    while stop == 0 && v > 0
        v = c(v);
        r = v==j;
        if r || any(q==v) %terminate if reachability or cycle is found
            stop = 1;
        end
        q = [q v];
    end