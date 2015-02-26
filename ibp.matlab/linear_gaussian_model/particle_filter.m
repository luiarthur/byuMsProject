function Zparticles = particle_filter(X,num_particles,true_alpha,true_sigma_x, true_sigma_A)

SAVEPARTIALRESULTS = 1;

Zparticles = cell(num_particles,1);

N = size(X,1);
time_1_obs = 0;
total_time = 0;
for n = 1:N
    particle_weight = ones(num_particles,1)/num_particles;
    if(n==1)
        disp(['Particle Filter:: Row ' num2str(n) '/' num2str(N) ]);
    else
        total_time = total_time + time_1_obs;
        disp(['Particle Filter:: Row ' num2str(n) '/' num2str(N) ', Elaps. Time: ' secs2hmsstr(total_time) ', Rem. Time: ' secs2hmsstr((total_time/n)*N-total_time)]);
    end
    tic

    if(n==1)
        for pind=1:num_particles
            num_first_dishes = poissrnd(true_alpha);
            Zparticles{pind} = ones(1,num_first_dishes);
            if (num_first_dishes==0)
                Zparticles{pind} = zeros(1);
            end
            Z_local = Zparticles{pind};
            particle_weight(pind) = logPX(X(n,:),Z_local(n,:),true_sigma_x, true_sigma_A, true_alpha);
        end
        logmax = max(particle_weight);
        particle_weight = exp(particle_weight-logmax);
        particle_weight = particle_weight/sum(particle_weight);
    else
        for pind=1:num_particles
            local_Z = Zparticles{pind};
            if(size(local_Z,1)==1)
                m_k = local_Z;
            else
                m_k = sum(local_Z);
            end
            shared_choices = m_k/n > rand(size(m_k));
            new_choices = poissrnd(true_alpha/n);

            Zparticles{pind} = [[Zparticles{pind} zeros(size(Zparticles{pind},1),new_choices)]; [shared_choices ones(1,new_choices)]];

            Z_local = Zparticles{pind};
            
            M = eye(n) - Z_local*inv(Z_local'*Z_local + (true_sigma_x^2 / true_sigma_A^2) * eye(size(Z_local,2)))*Z_local';
            
%             
%             S = M/true_sigma_x^2;
%             
%             ACCB = inv(S)
        
            
            A = eye(n)/true_sigma_x^2;
            Ainv = eye(n)*true_sigma_x^2;
            
            B = -inv(Z_local'*Z_local + (true_sigma_x^2 / true_sigma_A^2) * eye(size(Z_local,2)))/true_sigma_x^2;
            Binv = (-(Z_local'*Z_local + (true_sigma_x^2 / true_sigma_A^2) * eye(size(Z_local,2)))*true_sigma_x^2);
            Xli = Z_local;
            
            C_N1_inv = Ainv-Ainv*Xli*inv((Binv+Xli'*Ainv*Xli))*Xli'*Ainv;
            C_N1 = A+Xli*B*Xli';
            
            C_N_inv = C_N1(1:end-1,1:end-1)-1/C_N1(end,end)*C_N1(1:end-1,end)*C_N1(end,1:end-1);
            
             ACCB = C_N1_inv;
             C = ACCB(1:end-1,end);
             B = ACCB(end,end);
%             A = ACCB(1:end-1,1:end-1);
            x = X(1:n-1,:);
            y = X(n,:);
%              iA = inv(A)
%             mean = C'*iA*x;
%             covar = B - C'*iA*C;
             
             mean = C'*C_N_inv*x;
             covar = B - C'*C_N_inv*C;

            particle_weight(pind) = (-(y-mean)*inv(covar)*(y-mean)'/2)-0.5*length(y)*(log(2*pi*covar));%logPX(X(n,:),Z_local(n,:),true_sigma_x, true_sigma_A, true_alpha);
        end
        logmax = max(particle_weight);
        particle_weight = exp(particle_weight-logmax);

        particle_weight = particle_weight/sum(particle_weight);



    end

    Zparticles = resample(Zparticles,particle_weight,num_particles);

    time_1_obs = toc;
    if(SAVEPARTIALRESULTS==1 && mod(n,100)==0)
        fn = ['PF-out-' num2str(n)];
        save(fn)
    end
end



