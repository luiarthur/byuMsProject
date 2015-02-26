function [Z_sample, lP_sample, K_sample, alpha_sample, sigma_x_sample, sigma_A_sample] = hyper_sampler(X,num_samples,true_Z,true_alpha,true_sigma_x, true_sigma_A)


% build data structures to hold samples
Z_sample = cell(num_samples,1);
lP_sample = zeros(num_samples,1);
K_sample = zeros(num_samples,1); %this can always be computed as size(Z_sample{i},2)
sigma_x_sample = zeros(num_samples,1);
sigma_A_sample = zeros(num_samples,1);
alpha_sample = zeros(num_samples,1);

structure_error = zeros(num_samples,1);
in_degree_error = zeros(num_samples,1);

%initialize sampler with smart guesses
Z_sample{1} = round(rand(size(X,1),1));%size(true_Z)); % no links
sigma_x_sample(1) = .5;%true_sigma_x; %betarnd(1,1); % Beta(1,1) prior.  This could be rnd() as it Beta(1,1) is uniform over [0 1]
alpha_sample(1) = true_alpha;%gamrnd(1,1);   %<--- to ensure a Zsample in the beginning gamrnd(1,1); %Gamma(1,1) prior
sigma_A_sample(1) = .5;%true_sigma_A; %rand(); % uniform over [0 1] prior


K_sample(1) = size(Z_sample{1},2);

lP_sample(1) = logPXZ(X,Z_sample{1},sigma_x_sample(1),sigma_A_sample(1),alpha_sample(1));

% there are two metropolis (vs. Gibbs) samplings steps per sweep, one
% for lambda and one for epsilon -- make variables to keep track of
% the acceptance ratio so as to scale the proposal variance according to
% that ratio

num_sigma_x_moves_accepted = 0;
sigma_x_move_acceptance_ratio = zeros(num_samples,1);
sigma_x_move_acceptance_ratio(1) = 0;
sigma_x_proposal_variance = .01;
sigma_x_bound = [0 inf];

num_sigma_A_moves_accepted = 0;
sigma_A_move_acceptance_ratio = zeros(num_samples,1);
sigma_A_move_acceptance_ratio(1) = 0;
sigma_A_proposal_variance = .01;
sigma_A_bound = [0 inf];

num_alpha_moves_accepted = 0;
alpha_move_acceptance_ratio = zeros(num_samples,1);
alpha_move_acceptance_ratio(1) = 0;
alpha_proposal_variance = .1;
alpha_bound = [0 inf];

GRAPHICS = 0;

% if the truth is provided, compute the ground truth score
if(nargin > 2)
    true_model_score = logPXZ(X,true_Z,true_sigma_x,true_sigma_A,true_alpha);
end

Hn = sum(1./1:size(X,1));



for sweep = 2:num_samples 
    if(mod(sweep,2)==0)
        disp(['Sweep ' num2str(sweep) '/' num2str(num_samples) ]);
    end

    Z_sample{sweep} = sampZ(X,Z_sample{sweep-1},sigma_x_sample(sweep-1),sigma_A_sample(sweep-1),alpha_sample(sweep-1));
%     [Z_sample{sweep} ,Y_sample{sweep}] = cannonize(Z_sample{sweep},Y_sample{sweep});
%     [Z_sample{sweep} ,Y_sample{sweep}] = clean(Z_sample{sweep},Y_sample{sweep});



    % epsilon metropolis step
    sigma_A_proposal = sigma_A_sample(sweep-1)+randn*sigma_A_proposal_variance;
    if(sigma_A_proposal > sigma_A_bound(1) && sigma_A_proposal < sigma_A_bound(2))
        
        
        lp_sigma_A_proposal = logPXZ(X,Z_sample{sweep},sigma_x_sample(sweep-1),sigma_A_proposal,alpha_sample(sweep-1));
        lp_sigma_A = logPXZ(X,Z_sample{sweep},sigma_x_sample(sweep-1),sigma_A_sample(sweep-1),alpha_sample(sweep-1));

        log_acceptance_ratio = lp_sigma_A_proposal - lp_sigma_A;

        if(log(rand)<min(log_acceptance_ratio,0)) %accept the switch, otherwise stay pat..
            sigma_A_sample(sweep) = sigma_A_proposal;
            num_sigma_A_moves_accepted = num_sigma_A_moves_accepted+1;
        else
            sigma_A_sample(sweep) = sigma_A_sample(sweep-1);
        end
    else
        sigma_A_sample(sweep) = sigma_A_sample(sweep-1);
    end
    sigma_A_move_acceptance_ratio(sweep) = num_sigma_A_moves_accepted/sweep;
    %
    %     if(sweep > const1 & mod(sweep,5)==0)
    %         if(sigma_A_move_acceptance_ratio(sweep) < .2)
    %             sigma_A_proposal_variance = sigma_A_proposal_variance*0.9;
    %         elseif(sigma_A_move_acceptance_ratio(sweep) > .3)
    %             sigma_A_proposal_variance = sigma_A_proposal_variance*1.1;
    %         end
    %         if(sigma_A_proposal_variance>.5)
    %             sigma_A_proposal_variance=.5;
    %         end
    %     end


    % sigma_x metropolis step
    sigma_x_proposal = sigma_x_sample(sweep-1)+randn*sigma_x_proposal_variance;
    if(sigma_x_proposal > sigma_x_bound(1) && sigma_x_proposal < sigma_x_bound(2))
        lp_sigma_x_proposal = logPXZ(X,Z_sample{sweep},sigma_x_proposal,sigma_A_sample(sweep),alpha_sample(sweep-1));
        lp_sigma_x = logPXZ(X,Z_sample{sweep},sigma_x_sample(sweep-1),sigma_A_sample(sweep),alpha_sample(sweep-1));

        log_acceptance_ratio = lp_sigma_x_proposal - lp_sigma_x;

        if(log(rand)<min(log_acceptance_ratio,0)) %accept the switch, otherwise stay pat..
            sigma_x_sample(sweep) = sigma_x_proposal;
            num_sigma_x_moves_accepted = num_sigma_x_moves_accepted+1;
        else
            sigma_x_sample(sweep) = sigma_x_sample(sweep-1);
        end
    else
        sigma_x_sample(sweep) = sigma_x_sample(sweep-1);
    end
    %
    sigma_x_move_acceptance_ratio(sweep) = num_sigma_x_moves_accepted/sweep;
    %
    %     if(sweep > const1 & mod(sweep,5)==0)
    %         if(sigma_x_move_acceptance_ratio(sweep) < .2)
    %             sigma_x_proposal_variance = sigma_x_proposal_variance*0.9;
    %         elseif(sigma_x_move_acceptance_ratio(sweep) > .3)
    %             sigma_x_proposal_variance = sigma_x_proposal_variance*1.1;
    %         end
    %         if(sigma_x_proposal_variance>.5)
    %             sigma_x_proposal_variance=.5;
    %         end
    %     end


    % alpha gibbs step
%     K_plus = size(Z_sample{sweep},2);
%     alpha_sample(sweep) = gamrnd(1+K_plus,1/(1+Hn));
    
     % sigma_x metropolis step
    alpha_proposal = alpha_sample(sweep-1)+randn*alpha_proposal_variance;
    if(alpha_proposal > alpha_bound(1) && alpha_proposal < sigma_x_bound(2))
        lp_alpha_proposal = logPXZ(X,Z_sample{sweep},sigma_x_sample(sweep),sigma_A_sample(sweep),alpha_proposal);
        lp_alpha = logPXZ(X,Z_sample{sweep},sigma_x_sample(sweep),sigma_A_sample(sweep),alpha_sample(sweep-1));

        log_acceptance_ratio = lp_alpha_proposal - lp_alpha;

        if(log(rand)<min(log_acceptance_ratio,0)) %accept the switch, otherwise stay pat..
            alpha_sample(sweep) = alpha_proposal;
            num_alpha_moves_accepted = num_alpha_moves_accepted+1;
        else
            alpha_sample(sweep) = alpha_sample(sweep-1);
        end
    else
        alpha_sample(sweep) = alpha_sample(sweep-1);
    end
    %
    alpha_move_acceptance_ratio(sweep) = num_alpha_moves_accepted/sweep;
    

    lP_sample(sweep) = logPXZ(X,Z_sample{sweep},sigma_x_sample(sweep),sigma_A_sample(sweep),alpha_sample(sweep));
    K_sample(sweep) =  size(Z_sample{sweep},2);

    % debugging graphics for those interested
    if(nargin > 3)
        [Ek, EZZt, cur_in_degree_error, cur_structure_error] = inferstats(Z_sample(2:sweep),true_Z);
        structure_error(sweep) = cur_structure_error;
        in_degree_error(sweep) = cur_in_degree_error;
    else
        [Ek, EZZt] = inferstats(Z_sample(2:sweep));
    end



    if((mod(sweep,5)==0 || sweep == num_samples) && GRAPHICS ==1 )

        % sampler progress stuff
        figure(1)

        subplot(2,4,1)
        plot(lP_sample(2:sweep))
        if(nargin > 3)
            hold on
            plot(2:sweep,repmat(true_model_score,1,sweep-1),'--')
            hold off
        end
        title('Score')

        subplot(2,4,2)
        eh = plot(sigma_A_move_acceptance_ratio(2:sweep),'r');
        hold on
        el = plot(sigma_x_move_acceptance_ratio(2:sweep),'k');
        hold off
        legend([eh el],'\sigma_A','\sigma_x');
        title('Metropolis acceptance ratio')

        subplot(2,4,3)
        plot(2:sweep,sigma_A_sample(2:sweep))
        if(nargin>5)
            hold on
            plot(2:sweep,repmat(true_sigma_A,1,sweep-1),'--')
            hold off
        end
        title('\sigma_A')

        subplot(2,4,4)
        plot(2:sweep,sigma_x_sample(2:sweep))
        if(nargin>4)
            hold on
            plot(2:sweep,repmat(true_sigma_x,1,sweep-1),'--')
            hold off
        end
        title('\sigma_x')

        subplot(2,4,5)
        plot(2:sweep,alpha_sample(2:sweep))
        if(nargin>3)
            hold on
            plot(2:sweep,repmat(true_alpha,1,sweep-1),'--')
            hold off
        end
        title('\alpha')

%         subplot(2,4,6)
       

        subplot(2,4,7)
        plot(2:sweep,K_sample(2:sweep))
        if(nargin>3)
            hold on
            plot(2:sweep,repmat(size(true_Z,2),1,sweep-1),'--')
            hold off
        end
        title('K (current sample)')
        if(nargin>3)
            subplot(2,4,8)
            eh = plot(structure_error(2:sweep),'r');
            hold on
            el = plot(in_degree_error(2:sweep),'k');
            hold off
            legend([eh el],'Structure','In Degree');
            title('Error')
        end


        %posterior distribution visualizations
        figure(2)

        subplot(2,4,1)
        hist(sigma_A_sample(2:sweep),0:.05:1)
        title('\sigma_A')
        set(gca,'XLim',[0 1])

%         subplot(2,4,2)


        subplot(2,4,3)
        hist(sigma_x_sample(2:sweep),0:.05:1)
        title('\sigma_x')
        set(gca,'XLim',[0 1])

        subplot(2,4,5)
        hist(alpha_sample(2:sweep),0:.25:15)
        title('\alpha')
        set(gca,'XLim',[0 15])

        subplot(2,4,6)
        hist(K_sample(2:sweep),0:20)
        set(gca,'XLim',[0 20])
        if(nargin>3)
            title(sprintf('p(K|...), true K=%d ',size(true_Z,2)))
        else
            title('p(K|...)')
        end
        subplot(2,4,7)
        imagesc(EZZt);
        colormap(hot);
        title('E[ZZ'']')
        colorbar
        if(nargin>3)
            subplot(2,4,8)
            imagesc(true_Z*true_Z');
            colormap(hot);
            title('True ZZ''')
            colorbar
        end
        drawnow
    end

end