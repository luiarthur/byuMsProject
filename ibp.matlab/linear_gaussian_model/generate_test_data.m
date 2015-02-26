function [X,Z,A, sigma_X] = generate_test_data(num_sample_images)
if(nargin<1)
    num_sample_images = 1000;
end
image_part_1 = [0 1 0 0 0 0; 1 1 1 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
image_part_2 = [0 0 0 1 1 1; 0 0 0 1 0 1; 0 0 0 1 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
image_part_3 = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 1 0 0 0 0 0; 1 1 0 0 0 0; 1 1 1 0 0 0];
image_part_4 = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 1 1; 0 0 0 0 1 0; 0 0 0 0 1 0];
subplot(1,4,1)
imagesc(image_part_1)
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap('hot')
subplot(1,4,2)
imagesc(image_part_2)
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap('hot')
subplot(1,4,3)
imagesc(image_part_3)
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap('hot')
subplot(1,4,4)
imagesc(image_part_4)
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap('hot')

A = [reshape(image_part_1,1,numel(image_part_1)) ; ...
     reshape(image_part_2,1,numel(image_part_2)) ; ...
     reshape(image_part_3,1,numel(image_part_3)) ; ...
     reshape(image_part_4,1,numel(image_part_4))];
 
 
 num_latent_image_features = size(A,1);
 vector_image_size = size(A,2);
 sigma_X = .5;
 
Z = round(rand(num_sample_images,num_latent_image_features));

X = Z*A + randn(num_sample_images,vector_image_size)*sigma_X; % this should probably be sqrt(sigma_X) but the data looks different in that case 
% 
% figure
% for i=1:num_sample_images
%     subplot(round(sqrt(num_sample_images)), round(sqrt(num_sample_images)), i);
%     imagesc(reshape(X(i,:),size(image_part_1)));
%     colormap('hot');
% end

% now we have X and the true Z and A and true sigma_X