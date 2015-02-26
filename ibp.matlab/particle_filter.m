function [Zparticles, Yparticles] = particle_filter(X,Z,Y,alpha,epsilon,lambda,p,num_particles)
% function [Zparticles, Yparticles] =
% particle_filter(X,Z,Y,alpha,epsilon,lambda,p,num_particles)
% 
% X is the observed data matrix, Z and Y should be initialized to [].
% alpha, epsilon,lambda,p should be given reasonable values for the
% modelling domain.  num_particles controls how many particles are used in
% the SIS particle filter
%
%
GRAPHICS = 0;
Zparticles = cell(num_particles,1);
Yparticles = cell(num_particles,1);

N = size(X,1);
T = size(X,2);

X_one_inds = cell(N,1);
X_zero_inds = cell(N,1);

for(n=1:N)
    X_one_inds{n} = find(X(n,:)==1);
    X_zero_inds{n} = setdiff(1:T,X_one_inds{n});
end

%particle_weight = ones(num_particles,1)/num_particles;

total_time = 0;

for n = 1:N
    if(n==1)
        disp(['Particle Filter:: Row ' num2str(n) '/' num2str(N) ]);
    else
        total_time = total_time + time_1_obs;
        disp(['Particle Filter:: Row ' num2str(n) '/' num2str(N) ', Elaps. Time: ' secs2hmsstr(total_time) ', Rem. Time: ' secs2hmsstr((total_time/n)*N-total_time)]);
    end
    tic
    
    particle_weight = ones(num_particles,1)/num_particles;



    if(n==1)
        for(pind=1:num_particles)
            num_first_dishes = poissrnd(alpha);
            Zparticles{pind} = ones(1,num_first_dishes);
            if(num_first_dishes==0)
                Zparticles{pind} = zeros(1);
            end


            lpXiT = sum(X(1,:))*(log(1-(1-epsilon)*((1-lambda*p)^num_first_dishes)));
            lpXiT = lpXiT+sum(1-X(1,:))*(log((1-epsilon)*((1-lambda*p)^num_first_dishes)));
            particle_weight(pind) = lpXiT;

        end
        logmax = max(particle_weight);
        particle_weight = exp(particle_weight-logmax);

        particle_weight = particle_weight/sum(particle_weight);

    else
        for(pind=1:num_particles)
            local_Z = Zparticles{pind};
            K = size(local_Z,2);
            if(size(local_Z,1)==1)
                m_k = local_Z;
            else
                m_k = sum(local_Z);
            end
            shared_choices = m_k/n > rand(size(m_k));
            new_choices = poissrnd(alpha/n);

            Zparticles{pind} = [[Zparticles{pind} zeros(size(Zparticles{pind},1),new_choices)]; [shared_choices ones(1,new_choices)]];

            local_Y = Yparticles{pind};
            lpXiT = 0;
            %             for(i = 1:n-1)
            new_Z =  Zparticles{pind} ;

            e = new_Z(n,1:K)*local_Y;
            oneinds = X_one_inds{n};
            zeroinds = X_zero_inds{n};

            lpXiT = lpXiT+sum(log(1-(1-epsilon)*((1-lambda).^e(oneinds))*((1-lambda*p)^new_choices)));
            lpXiT = lpXiT+sum(log((1-epsilon)*((1-lambda).^e(zeroinds))*((1-lambda*p)^new_choices)));
                %     lpK1i = lpXiT -alpha/N + (K+newK)*log(alpha/N) - gammaln(K+newK+1);


            particle_weight(pind) = lpXiT;
        end
        logmax = max(particle_weight);
        particle_weight = exp(particle_weight-logmax);

        particle_weight = particle_weight/sum(particle_weight);



    end
%     figure(5)
% 
%         for(i=1:10)
%             for(j=1:10)
%                 subplot(10,10,(i-1)*10+j)
%                 imagesc(Zparticles{i*j})
%                 colormap('hot')
%             end
%         end
%         drawnow
%         pause(2)

    [Zparticles,Yparticles,indexes] = resample(Zparticles,Yparticles,particle_weight,num_particles);

%         for(i=1:10)
%             for(j=1:10)
%                 subplot(10,10,(i-1)*10+j)
%                 imagesc(Zparticles{i*j})
%                 colormap('hot')
%             end
%         end
%         drawnow
%         pause(.1)

    for(partind = 1:num_particles)
        num_new_rows = size(Zparticles{partind},2)-size(Yparticles{partind},1);
        loc_Z = Zparticles{partind};
        Yparticles{partind} = sampY_newrows_only(X(1:n,:),[ Yparticles{partind};zeros(num_new_rows,T) ],loc_Z,alpha,epsilon,lambda,p,size(Yparticles{partind},1)+1);
%         Yparticles{partind} = sampY(X(n,:),[ Yparticles{partind};zeros(num_new_rows,T) ],loc_Z(n,:),alpha,epsilon,lambda,p);
        %         Yparticles{partind} = sampY(X(1:n,:),[Yparticles{partind}; zeros(num_new_rows,T) ],Zparticles{partind},alpha,epsilon,lambda,p);
    end
    time_1_obs = toc;
end


% uinds = unique(indexes);
% num_unique_resampled_particles = length(uinds);
%
% for(urpi = 1:num_unique_resampled_particles)
%     shared_particle_indexes = find(indexes==uinds(urpi));
%     first_unique_particle_index = min(shared_particle_indexes);
%
%     Ysample =
%     sampYs(X,Zparticles{first_unique_particle_index},Yparticles{first_unique_particle_index},alpha,epsilon,lambda,p);
%
%     Yparticles{shared_particle_indexes


%     resample();
%
%     sample_Ys();



%     for(p=1:num_particles)
%
%         %     [Zsamples{i},Ysamphallucinate] = sampZ(X,Ysamples{i-1},Zsamples{i-1},alpha,epsilon,lambda,p,temp(current_temp_index));
%         %     Ysamples{i} = sampY(X,Ysamphallucinate,Zsamples{i},alpha,epsilon,lambda,p,temp(current_temp_index));
%         %     [Zsamples{i} ,Ysamples{i}] = cannonize(Zsamples{i},Ysamples{i});
%         %     [Zsamples{i} ,Ysamples{i}] = clean(Zsamples{i},Ysamples{i});
%         %     lP(i) = logPXYZ(X,Ysamples{i},Zsamples{i},alpha,epsilon,lambda,p,temp(current_temp_index));
%
%         for(t = 1:length(temp))
%             [Zsamptemp{t},Ysamptemphallucinate] = sampZ(X,Ysamptemp{t},Zsamptemp{t},alpha,epsilon,lambda,p,temp(t));
%             Ysamptemp{t} = sampY(X,Ysamptemphallucinate,Zsamptemp{t},alpha,epsilon,lambda,p,temp(t));
%             [Zsamptemp{t} ,Ysamptemp{t}] = cannonize(Zsamptemp{t},Ysamptemp{t});
%             [Zsamptemp{t} ,Ysamptemp{t}] = clean(Zsamptemp{t},Ysamptemp{t});
%             lPtemp(t) = logPXYZ(X,Ysamptemp{t},Zsamptemp{t},alpha,epsilon,lambda,p,temp(t));
%         end
%         Zsamples{i} = Zsamptemp{1};
%         Ysamples{i} = Ysamptemp{1};
%
%         lP(i) = lPtemp(1);
%         Kr(i) = size(Zsamples{i},2);
%
%         if(GRAPHICS)
%
%             figure(9)
%             subplot(2,1,1)
%             plot(lP(1:i))
%             title('log(P(X,Z,Y)) disregarding Y')
%             subplot(2,1,2)
%             plot(Kr(1:i))
%             title('K (latent dimensionality')
%             drawnow
%
%             figure(6)
%             imagesc(Zsamples{i})
%             title('Current Z Sample')
%
%             drawnow
%         end
%
%     end
end



function val = plusminus1rand()
val = round(rand)*2-1;
end