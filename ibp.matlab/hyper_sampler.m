function [Z_sample, Y_sample, lP_sample, K_sample, alpha_sample, epsilon_sample, lambda_sample, p_sample] = hyper_sampler(X,num_samples,true_Y,true_Z,true_alpha,true_epsilon,true_lambda,true_p)
% function [Z_sample, Y_sample, lP_sample, K_sample, alpha_sample,
% epsilon_sample, lambda_sample, p_sample] = hyper_sampler(X,num_samples,true_Y,true_Z,true_alpha,true_epsilon,true_lambda,true_p)
%
%  Generates 'num_samples' samples from the noisy-or IBP posterior for
%  given binary observation matrix X, all of the other parameters are used
%  only for diagnostic and optional plotting.  To enable diagnostic
%  plotting the function must be opened and GRAPHICS must be set to 1.
%
%  Z_sample and Y_sample are cell(num_samples,1) structures that contain
%  the samples.  lP_sample, K_sample, alpha_sample, epsilon_sample,
%  lambda_sample, p_sample are all [num_samples, 1] arrays containing
%  samples of the corresponding name 

GRAPHICS = 0;

% build data structures to hold samples
Z_sample = cell(num_samples,1);
Y_sample = cell(num_samples,1);
lP_sample = zeros(num_samples,1);
K_sample = zeros(num_samples,1); %this can always be computed as size(Z_sample{i},2)
p_sample = zeros(num_samples,1);
lambda_sample = zeros(num_samples,1);
epsilon_sample = zeros(num_samples,1);
alpha_sample = zeros(num_samples,1);

structure_error = zeros(num_samples,1);
in_degree_error = zeros(num_samples,1);

%initialize sampler with smart guesses
Z_sample{1} = zeros(size(X,1),1);%size(true_Z)); % no links
Y_sample{1} = zeros(size(Z_sample{1},2),size(X,2)); %with no links, this doesn't matter
p_sample(1) = betarnd(1,1); % Beta(1,1) prior.  This could be rnd() as it Beta(1,1) is uniform over [0 1]
alpha_sample(1) = gamrnd(1,1);   %<--- to ensure a Zsample in the beginning gamrnd(1,1); %Gamma(1,1) prior
lambda_sample(1) = rand(); % uniform over [0 1] prior
epsilon_sample(1) = rand();  %<--- to ensure a Zsample in the beginning rand(); %  ..

K_sample(1) = size(Z_sample{1},2);
lP_sample(1) = logPXYZ(X,Y_sample{1},Z_sample{1},alpha_sample(1),epsilon_sample(1),lambda_sample(1),p_sample(1));

% there are two metropolis (vs. Gibbs) samplings steps per sweep, one
% for lambda and one for epsilon -- make variables to keep track of
% the acceptance ratio so as to scale the proposal variance according to
% that ratio

num_lambda_moves_accepted = 0;
lambda_move_acceptance_ratio = zeros(num_samples,1);
lambda_move_acceptance_ratio(1) = 0;
lambda_proposal_variance = .05;
lambda_bound = [0 1];
num_epsilon_moves_accepted = 0;
epsilon_move_acceptance_ratio = zeros(num_samples,1);
epsilon_move_acceptance_ratio(1) = 0;
epsilon_proposal_variance = .05;
epsilon_bound = [0 1];

% if the truth is provided, compute the ground truth score
if(nargin > 2)
    true_model_score = logPXYZ(X,true_Y,true_Z,true_alpha,true_epsilon,true_lambda,true_p);
end

Hn = sum(1./([1:size(X,1)]));

const1 = 20;

for(sweep = 2:num_samples)
    if(mod(sweep,500)==0)
        disp(['Sweep ' num2str(sweep) '/' num2str(num_samples) ]);
    end


    [Z_sample{sweep},Y_sample{sweep}] = sampZ(X,Y_sample{sweep-1},Z_sample{sweep-1},alpha_sample(sweep-1),epsilon_sample(sweep-1),lambda_sample(sweep-1),p_sample(sweep-1));
    Y_sample{sweep} = sampY(X,Y_sample{sweep},Z_sample{sweep},alpha_sample(sweep-1),epsilon_sample(sweep-1),lambda_sample(sweep-1),p_sample(sweep-1));
    [Z_sample{sweep} ,Y_sample{sweep}] = cannonize(Z_sample{sweep},Y_sample{sweep});
    [Z_sample{sweep} ,Y_sample{sweep}] = clean(Z_sample{sweep},Y_sample{sweep});



    % epsilon metropolis step
    epsilon_proposal = epsilon_sample(sweep-1)+randn*epsilon_proposal_variance;
    if(epsilon_proposal > epsilon_bound(1) & epsilon_proposal < epsilon_bound(2))
        lp_epsilon_proposal = logPXYZ(X,Y_sample{sweep},Z_sample{sweep},alpha_sample(sweep-1),epsilon_proposal,lambda_sample(sweep-1),p_sample(sweep-1));
        lp_epsilon = logPXYZ(X,Y_sample{sweep},Z_sample{sweep},alpha_sample(sweep-1),epsilon_sample(sweep-1),lambda_sample(sweep-1),p_sample(sweep-1));

        log_acceptance_ratio = lp_epsilon_proposal - lp_epsilon;

        if(log(rand)<min(log_acceptance_ratio,0)) %accept the switch, otherwise stay pat..
            epsilon_sample(sweep) = epsilon_proposal;
            num_epsilon_moves_accepted = num_epsilon_moves_accepted+1;
        else
            epsilon_sample(sweep) = epsilon_sample(sweep-1);
        end
    else
        epsilon_sample(sweep) = epsilon_sample(sweep-1);
    end
    epsilon_move_acceptance_ratio(sweep) = num_epsilon_moves_accepted/sweep;
    %
    %     if(sweep > const1 & mod(sweep,5)==0)
    %         if(epsilon_move_acceptance_ratio(sweep) < .2)
    %             epsilon_proposal_variance = epsilon_proposal_variance*0.9;
    %         elseif(epsilon_move_acceptance_ratio(sweep) > .3)
    %             epsilon_proposal_variance = epsilon_proposal_variance*1.1;
    %         end
    %         if(epsilon_proposal_variance>.5)
    %             epsilon_proposal_variance=.5;
    %         end
    %     end


    % lambda metropolis step
    lambda_proposal = lambda_sample(sweep-1)+randn*lambda_proposal_variance;
    if(lambda_proposal > lambda_bound(1) & lambda_proposal < lambda_bound(2))
        lp_lambda_proposal = logPXYZ(X,Y_sample{sweep},Z_sample{sweep},alpha_sample(sweep-1),epsilon_sample(sweep),lambda_proposal,p_sample(sweep-1));
        lp_lambda = logPXYZ(X,Y_sample{sweep},Z_sample{sweep},alpha_sample(sweep-1),epsilon_sample(sweep),lambda_sample(sweep-1),p_sample(sweep-1));

        log_acceptance_ratio = lp_lambda_proposal - lp_lambda;

        if(log(rand)<min(log_acceptance_ratio,0)) %accept the switch, otherwise stay pat..
            lambda_sample(sweep) = lambda_proposal;
            num_lambda_moves_accepted = num_lambda_moves_accepted+1;
        else
            lambda_sample(sweep) = lambda_sample(sweep-1);
        end
    else
        lambda_sample(sweep) = lambda_sample(sweep-1);
    end
    %
    lambda_move_acceptance_ratio(sweep) = num_lambda_moves_accepted/sweep;
    %
    %     if(sweep > const1 & mod(sweep,5)==0)
    %         if(lambda_move_acceptance_ratio(sweep) < .2)
    %             lambda_proposal_variance = lambda_proposal_variance*0.9;
    %         elseif(lambda_move_acceptance_ratio(sweep) > .3)
    %             lambda_proposal_variance = lambda_proposal_variance*1.1;
    %         end
    %         if(lambda_proposal_variance>.5)
    %             lambda_proposal_variance=.5;
    %         end
    %     end




    % p gibbs step
    num_ones_in_Y = sum(sum(Y_sample{sweep}));
    num_zeros_in_Y = prod(size(Y_sample{sweep})) - num_ones_in_Y;

    p_sample(sweep) = betarnd(num_ones_in_Y+1,num_zeros_in_Y+1);

    % alpha gibbs step
    K_plus = size(Z_sample{sweep},2);
    alpha_sample(sweep) = gamrnd(1+K_plus,1/(1+Hn));


    lP_sample(sweep) = logPXYZ(X,Y_sample{sweep},Z_sample{sweep},alpha_sample(sweep),epsilon_sample(sweep),lambda_sample(sweep),p_sample(sweep));
    K_sample(sweep) =  size(Z_sample{sweep},2);

    % debugging graphics for those interested
    if(nargin > 3)
        [Ek, EZZt, cur_in_degree_error, cur_structure_error] = inferstats(Z_sample(2:sweep),true_Z);
        structure_error(sweep) = cur_structure_error;
        in_degree_error(sweep) = cur_in_degree_error;
    else
        [Ek, EZZt] = inferstats(Z_sample(2:sweep));
    end



    if((mod(sweep,250)==0 || sweep == num_samples) && GRAPHICS)

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
        eh = plot(epsilon_move_acceptance_ratio(2:sweep),'r');
        hold on
        el = plot(lambda_move_acceptance_ratio(2:sweep),'k');
        hold off
        legend([eh el],'\epsilon','\lambda');
        title('Metropolis acceptance ratio')

        subplot(2,4,3)
        plot(2:sweep,epsilon_sample(2:sweep))
        if(nargin>5)
            hold on
            plot(2:sweep,repmat(true_epsilon,1,sweep-1),'--')
            hold off
        end
        title('\epsilon')

        subplot(2,4,4)
        plot(2:sweep,lambda_sample(2:sweep))
        if(nargin>6)
            hold on
            plot(2:sweep,repmat(true_lambda,1,sweep-1),'--')
            hold off
        end
        title('\lambda')

        subplot(2,4,5)
        plot(2:sweep,alpha_sample(2:sweep))
        if(nargin>4)
            hold on
            plot(2:sweep,repmat(true_alpha,1,sweep-1),'--')
            hold off
        end
        title('\alpha')

        subplot(2,4,6)
        plot(2:sweep,p_sample(2:sweep))
        if(nargin>7)
            hold on
            plot(2:sweep,repmat(true_p,1,sweep-1),'--')
            hold off
        end
        title('p')

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
        hist(lambda_sample(2:sweep),0:.05:1)
        title('\lambda')
        set(gca,'XLim',[0 1])

        subplot(2,4,2)
        hist(epsilon_sample(2:sweep),0:.05:1)
        title('\epsilon')
        set(gca,'XLim',[0 1])

        subplot(2,4,3)
        hist(p_sample(2:sweep),0:.05:1)
        title('p')
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