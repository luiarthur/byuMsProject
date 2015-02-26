GRAPHICS = 1;
num_samples = 250;
Zsamples = round(rand([size(Z) num_samples]));
Ysamples = round(rand([size(Y) num_samples]));
% Zsamples(:,:,1) = Z;
%Ysamples(:,:,1) = Y;
lP = zeros(num_samples(1),1);

lP(1) = logPXYZ(X,Ysamples(:,:,1),Zsamples(:,:,1),alpha,epsilon,lambda,p);

for(i = 2:num_samples)
    disp(['Iter ' num2str(i) '/' num2str(num_samples) ]);
    Zsamples(:,:,i) = sampZ(X,Ysamples(:,:,i-1),Zsamples(:,:,i-1),alpha,epsilon,lambda,p);
    %     Zsamples(:,:,i) = Zsamples(:,:,i-1);
    Ysamples(:,:,i) = sampY(X,Ysamples(:,:,i-1),Zsamples(:,:,i),alpha,epsilon,lambda,p);
%     Ysamples(:,:,i) = Ysamples(:,:,i-1);
    lP(i) = logPXYZ(X,Ysamples(:,:,i),Zsamples(:,:,i),alpha,epsilon,lambda,p);
    [Zsamples(:,:,i) ,Ysamples(:,:,i)] = cannonize(Zsamples(:,:,i),Ysamples(:,:,i));
    if(GRAPHICS)

        figure(9)
        plot(lP(1:i))



        pnX = 1-(1-lambda).^(Zsamples(:,:,i)*Ysamples(:,:,i))*(1-epsilon);
        flips = rand(size(pnX));
        nX=zeros(size(pnX));
        nX(find(flips<pnX))=1;

        figure(5)
        imagesc(nX)
        title(['X Sampled  log(P(X,Y,Z)) = ' num2str(lP(i))])

        figure(6)
        imagesc(Zsamples(:,:,i))
        title('Z Learned')
        figure(7)
        imagesc(Ysamples(:,:,i))
        title('Y Learned')

        drawnow
    end

end



pnX = 1-(1-lambda).^(Zsamples(:,:,end)*Ysamples(:,:,end))*(1-epsilon);
flips = rand(size(pnX));
nX=zeros(size(pnX));
nX(find(flips<pnX))=1;

figure(5)
imagesc(nX)
title(['X Sampled  log(P(X,Y,Z)) = ' num2str(lP(end))])
figure(6)
imagesc(Zsamples(:,:,end))
title('Z Learned')
figure(7)
imagesc(Ysamples(:,:,end))
title('Y Learned')