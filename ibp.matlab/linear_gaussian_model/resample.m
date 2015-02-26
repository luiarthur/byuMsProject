% Copyright October, 2005, Brown University, Providence, RI. 
% All Rights Reserved 

% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a commercial
% product is hereby granted without fee, provided that the above copyright
% notice appear in all copies and that both that copyright notice and this
% permission notice appear in supporting documentation, and that the name of
% Brown University not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission. 

% BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
% INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
% PARTICULAR PURPOSE. IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR ANY
% SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
% RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
% CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
% CONNECTION WITH THE USE.

% Author: Frank Wood fwood@cs.brown.edu.

function [retZparticles,indexes] = resample(Zparticles,w,N)    

retZparticles = Zparticles;
indexes = zeros(size(Zparticles));
cdf = cumsum(w);
p = rand(1,N);
p = sort(p);

picked = zeros(size(Zparticles));
for i=1:N 
    
    pind = min(find(cdf>p(i)));
    picked(pind) = picked(pind)+1;
    
%     if(cdf(cdfind)>p(i))
%         picked(cdfind) = picked(cdfind)+1;
%     else
%         
%         while(cdf(cdfind)<p(i))
%             cdfind=cdfind+1;
%         end
%         picked(cdfind) = picked(cdfind)+1;
%     end
    
end

rxind=1;
for i=1:N 
    if(picked(i)>0)
        for j=1:picked(i)
            retZparticles{rxind} = Zparticles{i};
            indexes(rxind) = i;
            rxind=rxind+1;
        end
    end
end

