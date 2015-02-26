ksum = 0;
for(i=1:num_samples)
   ksum = ksum+size(Zsamples{i},2); 
end
ksum/num_samples