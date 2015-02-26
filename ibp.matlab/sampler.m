function [Zsamples, Ysamples, lP, Kr, Tsamples] = sampler(X,Z,Y,alpha,epsilon,lambda,p,num_samples, use_simulated_tempering, start_from_truth)

if(nargin<10)
    start_from_truth = 0;
end
if(nargin<9)
    use_simulated_tempering = 0;
end

debug = 0;
GRAPHICS = 0;
Zsamples = cell(num_samples,1);
Ysamples = cell(num_samples,1);
Tsamples = zeros(num_samples,1);

ground_truth_score = logPXYZ(X,Y,Z,alpha,epsilon,lambda,p,1);

if(start_from_truth==1)
    Zsamples{1} = Z;
    Ysamples{1} = Y;
else
    Zsamples{1} = zeros(size(Z));

    Ysamples{1} = zeros(size(Z,2),size(X,2));%round(rand(size(Y)));
%     Ysamp = sampY(X,Ysamples{1},Zsamples{1},alpha,epsilon,lambda,p);
%     [Zsamples{1} ,Ysamples{1}] = cannonize(Zsamples{1},Ysamp);
end

lP = zeros(num_samples(1),1);
Kr = zeros(num_samples(1),1);



% temp = [1 2 3 4 5];

% temp = [1 .7071 .5 .25 .125];
if(use_simulated_tempering)
    temp(1) = 1/20;
else
    temp = 1;
end
% temp_prior = (temp)/norm( (temp));
% 
% for(t=1:length(temp))
%     temp_prior(t) = logPXYZ(X,Ysamples{1},Zsamples{1},alpha,epsilon,lambda,p,temp(t));
% end
% temp_prior = (1./exp(temp_prior/1000))/norm(1./exp(temp_prior/1000));
% log_temp_prior = log(temp_prior);
num_moves_attempted_ij = zeros(length(temp),length(temp));
num_moves_accepted_ij = zeros(size(num_moves_attempted_ij));

if(use_simulated_tempering)
    current_temp_index = length(temp);
else
    current_temp_index = 1;
end

Tsamples(1) = temp(1);
Kr(1) = size(Zsamples{1},2);
lP(1) = logPXYZ(X,Ysamples{1},Zsamples{1},alpha,epsilon,lambda,p,Tsamples(1) );

Zsamptemp = cell(1,length(temp));
Ysamptemp = cell(1,length(temp));
lPtemp = zeros(1,length(temp));
for(t=1:length(temp))
    Zsamptemp{t} = Zsamples{1};
    Ysamptemp{t} = Ysamples{1};
end
consider_even_swaps=1;
for(i = 1:num_samples)
    if(mod(i,500)==0)
        disp(['Iter ' num2str(i) '/' num2str(num_samples) ]);
    end
%     [Zsamples{i},Ysamphallucinate] = sampZ(X,Ysamples{i-1},Zsamples{i-1},alpha,epsilon,lambda,p,temp(current_temp_index));
%     Ysamples{i} = sampY(X,Ysamphallucinate,Zsamples{i},alpha,epsilon,lambda,p,temp(current_temp_index));
%     [Zsamples{i} ,Ysamples{i}] = cannonize(Zsamples{i},Ysamples{i});
%     [Zsamples{i} ,Ysamples{i}] = clean(Zsamples{i},Ysamples{i});
%     lP(i) = logPXYZ(X,Ysamples{i},Zsamples{i},alpha,epsilon,lambda,p,temp(current_temp_index));

    for(t = 1:length(temp))
        [Zsamptemp{t},Ysamptemphallucinate] = sampZ(X,Ysamptemp{t},Zsamptemp{t},alpha,epsilon,lambda,p,temp(t));
        Ysamptemp{t} = sampY(X,Ysamptemphallucinate,Zsamptemp{t},alpha,epsilon,lambda,p,temp(t));
        [Zsamptemp{t} ,Ysamptemp{t}] = cannonize(Zsamptemp{t},Ysamptemp{t});
        [Zsamptemp{t} ,Ysamptemp{t}] = clean(Zsamptemp{t},Ysamptemp{t});
        lPtemp(t) = logPXYZ(X,Ysamptemp{t},Zsamptemp{t},alpha,epsilon,lambda,p,temp(t));
    end
%     lPtemp
%     swaps = nchoosek(1:length(temp),2);
%     
%     if(consider_even_swaps==2)
%         consider_even_swaps = 1;
%     else
%         consider_even_swaps = 2;
%     end
%     
%     for(s=consider_even_swaps:2:size(swaps,1))
%         ii = swaps(s,1);
%         jj = swaps(s,2);
%         
%         hixj = logPXYZ(X,Ysamptemp{jj},Zsamptemp{jj},alpha,epsilon,lambda,p,temp(ii));
%         hjxi = logPXYZ(X,Ysamptemp{ii},Zsamptemp{ii},alpha,epsilon,lambda,p,temp(jj));
%         hixi = lPtemp(ii);
%         hjxj = lPtemp(jj);
%         
%         lr = (hixj+hjxi)-(hixi+hjxj);
%         disp(sprintf('Swap considered between temp indexes''s %d <-> %d: log(r) = %0.3f',ii,jj,lr))
% %         ar = 1/(1+exp((hixi+hjxj)-(hixj+hjxi)));
%         if(log(rand)<min(lr,0)) %swap
% %         if(rand <ar)
%             disp(sprintf('Swap accepted between temp indexes''s %d <-> %d',ii,jj))
%             z_temp = Zsamptemp{ii};
%             Zsamptemp{ii} = Zsamptemp{jj};
%             Zsamptemp{jj} = z_temp;
%             
%             y_temp = Ysamptemp{ii};
%             Ysamptemp{ii} = Ysamptemp{jj};
%             Ysamptemp{jj} = y_temp;
%             
%             lp_temp = lPtemp(ii);
%             lPtemp(ii) = lPtemp(jj);
%             lPtemp(jj) = lp_temp;
%             
%         end
%     end
    
    Zsamples{i} = Zsamptemp{1};
    Ysamples{i} = Ysamptemp{1};
    Tsamples(i) = temp(1);
    lP(i) = lPtemp(1);

    if(use_simulated_tempering)
    temp = 1/((1/temp)/(1.000001^i));
    if(temp>1)
        temp =1;
    end
    end

    
%     if(use_simulated_tempering)
%         Tsamples(i) = current_temp_index;
%         
% %         cti = hist(i,(1:length(temp))*5);
% %         current_temp_index = find(cti(end:-1:1)==1);
%         if(current_temp_index == length(temp))
%             proposed_next_temp_index = current_temp_index-1;
%             qij = 1;
%         elseif(current_temp_index == 1)
%             proposed_next_temp_index = 2;
%             qij = 1;
%         else
%             proposed_next_temp_index = current_temp_index+plusminus1rand();
%             qij = .5;
%         end
% 
%         if(proposed_next_temp_index== length(temp) | proposed_next_temp_index == 1)
%             qji = 1;
%         else
%             qji = .5;
%         end
% 
% 
% 
%         hix = lP(i);
%         hjx = logPXYZ(X,Ysamples{i},Zsamples{i},alpha,epsilon,lambda,p,temp(proposed_next_temp_index));
% 
%         lr = hjx + log_temp_prior(proposed_next_temp_index)+log(qji) - (hix + log_temp_prior(current_temp_index)+log(qij))
%         %         r = (exp(hjx)*temp_prior(proposed_next_temp_index)*qji) / (exp(hix)*temp_prior(current_temp_index)*qij)
% 
%         num_moves_attempted_ij(current_temp_index,proposed_next_temp_index) = num_moves_attempted_ij(current_temp_index,proposed_next_temp_index)+1;
% 
%         pre_move_temp_index = current_temp_index;
%         if(log(rand)<min(lr,0))
%             current_temp_index = proposed_next_temp_index;
%             num_moves_accepted_ij(pre_move_temp_index,proposed_next_temp_index) =num_moves_accepted_ij(pre_move_temp_index,proposed_next_temp_index)+1;
%         end
% 
%         acceptance_ratio = num_moves_accepted_ij(pre_move_temp_index,proposed_next_temp_index)/num_moves_attempted_ij(pre_move_temp_index,proposed_next_temp_index);
% 
% %         c0 = 500;
% %         n0 = 1;
% %         log_temp_prior(current_temp_index) = log_temp_prior(current_temp_index) -c0/(i/4+n0);
% %         otherinds = setdiff(1:length(temp_prior),current_temp_index);
% %         log_temp_prior(otherinds) = log_temp_prior(otherinds) + c0/(length(log_temp_prior)*(i/4+n0))
% %         re_normalized = exp(log_temp_prior)./sum(exp(log_temp_prior));
% %         re_normalized(find(re_normalized ==0)) = realmin;
% %         log_temp_prior = log(re_normalized);
% %         if(mod(i,20)==0)
% %             disp('foo');
% %         end
% 
%         %
%         %         if(acceptance_ratio < .4)
%         %             log_temp_prior(pre_move_temp_index) = log_temp_prior(pre_move_temp_index)+lr/2;
%         %             otherinds = setdiff(1:length(temp_prior),pre_move_temp_index);
%         %             log_temp_prior(otherinds) =  log_temp_prior(otherinds)-lr/(2*length(otherinds));
%         %             re_normalized = exp(log_temp_prior)./sum(exp(log_temp_prior));
%         %             re_normalized(find(re_normalized ==0)) = log(realmin);
%         %             log_temp_prior = log(re_normalized);
%         %         elseif(acceptance_ratio > .6)
%         %             log_temp_prior(proposed_next_temp_index) = log_temp_prior(proposed_next_temp_index)-10;
%         %             log_temp_prior(setdiff(1:length(temp_prior),current_temp_index)) =  log_temp_prior(setdiff(1:length(temp_prior),proposed_next_temp_index))+10;
%         % %             temp_prior = temp_prior./norm(temp_prior);
%         %             re_normalized = exp(log_temp_prior)./sum(exp(log_temp_prior));
%         %             re_normalized(find(re_normalized ==0)) = log(realmin);
%         %             log_temp_prior = log(re_normalized);
%         %         end
% 
%         %         if(rand<min(r,1))
%         %             current_temp_index = proposed_next_temp_index;
%         %         end
% 
%     end
% 

    Kr(i) = size(Zsamples{i},2);
    if(0) %GRAPHICS | 
        if(~use_simulated_tempering)
            figure(9)
            subplot(2,1,1)
            plot(lP(1:i))
            hold on
            plot(1:i,repmat(ground_truth_score,1,i),'--')
            hold off
            title('log(P(X,Z,Y)) disregarding Y')
            subplot(2,1,2)
            plot(Kr(1:i))
            title('K (latent dimensionality')
            drawnow
        else
            figure(9)
            subplot(3,1,1)
            plot(lP(1:i))
            title('log(P(X,Z,Y)) disregarding Y')
            subplot(3,1,2)
            plot(Kr(1:i))
            title('K (latent dimensionality')
            subplot(3,1,3)
            plot(1./Tsamples(1:i))
            title('Temperature Index (1 is cool)')
            drawnow
        end



        %         pnX = 1-(1-lambda).^(Zsamples{i}*Ysamples{i})*(1-epsilon);
        %         flips = rand(size(pnX));
        %         nX=zeros(size(pnX));
        %         nX(find(flips<pnX))=1;
        %
        %         figure(5)
        %         imagesc(nX)
        %         title(['X Sampled  log(P(X,Y,Z)) = ' num2str(lP(i))])

        figure(6)
        MCMCMC=0;
        if(MCMCMC)
            for(t=1:length(temp))
                subplot(length(temp),1,t);
                imagesc(Zsamptemp{t})
                title(sprintf('Temp index = %d',t))
            end
        else
        imagesc(Zsamples{i})
        title('Current Z Sample')
        end
        %         figure(7)
        %         imagesc(Ysamples{i})
        %         title('Current Y Sample')

        drawnow
    end

end

% [Ek, EZZt] = inferstats(Zsamples,0);
%
% if(GRAPHICS)
%
%     pnX = 1-(1-lambda).^(Zsamples{end}*Ysamples{end})*(1-epsilon);
%     flips = rand(size(pnX));
%     nX=zeros(size(pnX));
%     nX(find(flips<pnX))=1;
%
%     figure(5)
%     imagesc(nX)
%     title(['X Sampled  log(P(X,Y,Z)) = ' num2str(lP(end))])
%     figure(6)
%     imagesc(Zsamples{end})
%     title('Z Learned')
%     figure(7)
%     imagesc(Ysamples{end})
%     title('Y Learned')
% end

function val = plusminus1rand()
val = round(rand)*2-1;