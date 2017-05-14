% This Matlab program demomstrates the level set algorithm in paper:
%    "Active contours driven by local image fitting energy" 
%    to appear in Pattern Recognition, 2010
% Author: Kaihua Zhang, Huihui Song and Lei Zhang
% E-mail: zhkhua@mail.ustc.edu.cn, cslzhang@comp.polyu.edu.hk  
% http://www4.comp.polyu.edu.hk/~cslzhang/

%  Notes:
%   1. Some parameters may need to be modified for different types of images. Please contact the author if any problem regarding the choice of parameters.
%   2. Intial contour should be set properly.

% Date: 5/11/2009

c0 =2;
Img = imread('vessel2.bmp');
Img = double(Img(:,:,1));
I=Img(:,:,1);
[nrow,ncol] = size(Img);

   %Draw mask
       figure;imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;axis equal;
      % text(6,6,'Left click to get points, right click to get end point','FontSize',[12],'Color', 'g');
       BW=roipoly;
       phi=c0*2*(0.5-BW);
       hold on;
       [c,h]=contour(phi,[0 0],'r');
      set(h, 'linewidth', 2.5);
       hold off;
       %end draw map;

pause(0.01);
figure;  imagesc(Img,[0 255]);colormap(gray);hold on;
contour(phi,[0 0],'b');

sigma =3;% the key parameter which needs to be tuned properly.
sigma_phi = 0.5;% the variance of regularized Gaussian kernel
K = fspecial('gaussian',2*round(2*sigma)+1,sigma);
K_phi = fspecial('gaussian',5,sigma_phi);


timestep = 0.1;
epsilon = 0.5;
for n = 1:2000   
      [phi,f1,f2,Hphi]= LIF_2D(Img,phi,timestep,epsilon,K);
      phi = conv2(phi,K_phi,'same');
      if mod(n,40)==0
      pause(0.0001);
      imagesc(Img,[0 255]);colormap(gray)
      hold on;contour(phi,[0 0],'b');
      iterNum=[num2str(n), ' iterations'];        
      title(iterNum);        
      hold off;
      end
end

figure;
mesh(phi);