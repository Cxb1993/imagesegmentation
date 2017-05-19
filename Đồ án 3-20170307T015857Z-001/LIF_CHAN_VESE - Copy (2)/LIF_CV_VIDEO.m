function   LIV_CV_VDIEO(sigma_phi,num_inter,sigma,file_name,mu,a_1,timestep,epsilon)

workingDir = 'C:\Users\doantuananh\Downloads\LIF_CHAN_VESE - Copy (2)';
%mkdir(workingDir)
mkdir(workingDir,'output');
outputVideo = VideoWriter(fullfile(workingDir,'shuttle_output.avi'));
outputVideo.FrameRate = 15;

open(outputVideo)
 %% ?o?n này thêm
mov=VideoReader('shuttle_out2.avi');
%vidFrames=read(mov);
nFrames=mov.NumberOfFrames;
%vidFramesgrey(:,:,:)=0.2989*vidFrames(:,:,1,:)+0.5870*vidFrames(:,:,2,:)+0.1140*vidFrames(:,:,3,:);

%vidFramesgrey=double(vidFramesgrey);
% for i=1:1
%    imshow(vidFramesgrey(:,:,i),[]);  %frames are grayscale
% end
 %%?o?n này thêm
c0 =2;
%Img = imread(file_name);
%Img = double(Img(:,:,1));
sigma =3;% the key parameter which needs to be tuned properly.
sigma_phi = 0.45;% the variance of regularized Gaussian kernel
K = fspecial('gaussian',2*round(2*sigma)+1,sigma);
K_phi = fspecial('gaussian',5,sigma_phi);
%I=Img(:,:,1);
%[nrow,ncol] = size(Img);

%    %Draw mask
%         figure;imagesc(vidFramesgrey(:,:,1), [0, 255]);colormap(gray);hold on; axis off;axis equal;
%        % text(6,6,'Left click to get points, right click to get end point','FontSize',[12],'Color', 'g');
%         BW=roipoly;
%         phi=c0*2*(0.5-BW);
%         hold on;
%         [c,h]=contour(phi,[0 0],'r');
%        set(h, 'linewidth', 2.5);
%         hold off;
%        %end draw map;
% 
%  pause(0.01);
%  figure;  imagesc(vidFramesgrey(:,:,1),[0 255]);colormap(gray);hold on;
%  contour(phi,[0 0],'b');
% 
timestep = 0.05;
epsilon = 1.0;
% % for m=2:40
% %     for n = 1:num_inter
% %         old=phi;
% %         [force,phi,f1,f2,Hphi]= Caculate_next_phi(mu,a_1,double(vidFramesgrey(:,:,m)),phi,timestep,epsilon,K);
% %         
% %         phi = conv2(phi,K_phi,'same');
% %         new=phi;
% %         indicator=checkstop(old,new,timestep);
% %                if(indicator)
% %         %           imagesc(Img,[0 255]);colormap(gray)
% %         %       hold on;[c,h]=contour(phi,[0 0],'b');
% %         %       set(h, 'linewidth', 2.5);
% %         %       iterNum=[num2str(n), ' iterations'];
% %         %       title(iterNum);
% %         %       hold off;
% %         %       figure;
% %         %       mesh(phi);
% %         %       return;
% %             pause(0.0001);
% %             imagesc(vidFramesgrey(:,:,m),[0 255]);colormap(gray)
% %             hold on;[c,h]=contour(phi,[0 0],'b');
% %             set(h, 'linewidth', 2.5);
% %             iterNum=[num2str(n), ' iterations'];
% %             title(iterNum);
% %             hold off;
% %                break;
% %          end;
% %         %if mod(n,num_inter)==0
% %             pause(0.0001);
% %             imagesc(vidFramesgrey(:,:,m),[0 255]);colormap(gray)
% %             hold on;[c,h]=contour(phi,[0 0],'b');
% %             set(h, 'linewidth', 2.5);
% %             iterNum=[num2str(n), ' iterations'];
% %             title(iterNum);
% %             hold off;
% %         %end
% %         
% %     end
% %     img1=getframe(gcf);
% %     writeVideo(outputVideo,img1);
% % end
% % figure;
% % mesh(phi);
%viet output


%close(outputVideo)
phi=0;
videofig(mov.NumberOfFrames, @(frm) redraw(frm, mov,mu,a_1,phi,timestep,epsilon,K,num_inter,K_phi,c0));

% Display initial frame
redraw(1, mov,mu,a_1,phi,timestep,epsilon,K,num_inter,K_phi,c0);
end
function redraw(frame, vidObj,mu,a_1,phi,timestep,epsilon,K,num_inter,K_phi,c0)
% REDRAW  Process a particular frame of the video
%   REDRAW(FRAME, VIDOBJ)
%       frame  - frame number to process
%       vidObj - VideoReader object

% Read frame
f = vidObj.read(frame);

        
vidFramesgrey(:,:)=0.2989*f(:,:,1)+0.5870*f(:,:,2)+0.1140*f(:,:,3);

if frame==1
    figure;imagesc(vidFramesgrey(:,:), [0, 255]);colormap(gray);hold on; axis off;axis equal;
       % text(6,6,'Left click to get points, right click to get end point','FontSize',[12],'Color', 'g');
        BW=roipoly;
        phi0=c0*2*(0.5-BW);
        hold on;
        [c,h]=contour(phi0,[0 0],'r');
       set(h, 'linewidth', 2.5);
        hold off;
        phi=phi0;
        for n = 1:num_inter
        old=phi;
        [force,phi,f1,f2,Hphi]= Caculate_next_phi(mu,a_1,double(vidFramesgrey(:,:)),phi,timestep,epsilon,K);
         phi = conv2(phi,K_phi,'same');
        new=phi;
        indicator=checkstop(old,new,timestep);
               if(indicator)
                   break;
               end
        end
else
     for n = 1:num_inter
        old=phi;
        [force,phi,f1,f2,Hphi]= Caculate_next_phi(mu,a_1,double(vidFramesgrey(:,:)),phi,timestep,epsilon,K);
         phi = conv2(phi,K_phi,'same');
        new=phi;
        indicator=checkstop(old,new,timestep);
               if(indicator)
                   break;
               end
        end
end

% Overlay edge on original image
%f3 = bsxfun(@plus, f,  phi);

% Display
image(vidFramesgrey);
hold on;
 [c,h]=contour(phi,[0 0],'r');
 hold off;
axis image off

end


function [force,phi,f1,f2,Hphi] = Caculate_next_phi(mu,a_1,I,phi,timestep,epsilon,K)
inidx = find(phi>=0); % frontground index
outidx = find(phi<0); % background index
phi = NeumannBoundCond(phi);
P = double(I);
L = im2double(P(:,:,1)); % get one image component
force_image = 0; %
c1 = sum(sum(L.*Heaviside1(phi)))/(length(inidx)+eps); % average inside of Phi0
c2 = sum(sum(L.*(1-Heaviside1(phi))))/(length(outidx)+eps); % verage outside of Phi0

force_image=-(L-c1).^2+(L-c2).^2+force_image; 
 force = force_image;
 %force = force_image;
 %force = force./(max(max(abs(force))));
Hphi = Heaviside(phi,epsilon);
DiracPhi = Delta(phi,epsilon);
%penalizeTerm=mu*(4*del2(phi)-kappa(phi));
[f1,f2] = Local_Avr(I,Hphi,K);

phi = phi + timestep*(DiracPhi.*((a_1.*((I - f1.*Hphi - f2.*(1 - Hphi)).*(f1 - f2)))+(1-a_1)*force));
end
function H = Heaviside(phi,epsilon)
H = 0.5*(1+(2/pi)*atan(phi./epsilon));
end
function Delta_h = Delta(phi,epsilon)
Delta_h = (epsilon/pi)./(epsilon^2+phi.^2);
end

function [f1,f2] = Local_Avr(I,H,K)

f1 = conv2(I.*H,K,'same');
c1 = conv2(H,K,'same');
f1 = f1./c1;
f2 = conv2(I.*(1-H),K,'same');
c2 = conv2(1-H,K,'same');
f2 = f2./c2;
end
function k = curvature_central(u)
% compute curvature for u with central difference scheme
[ux,uy] = gradient(u);
normDu = sqrt(ux.^2+uy.^2+1e-10);
Nx = ux./normDu;
Ny = uy./normDu;
[nxx,junk] = gradient(Nx);
[junk,nyy] = gradient(Ny);
k = nxx+nyy;
end
function H=Heaviside1(z)
% Heaviside step function (smoothed version)
% Copyright (c) 2009, 
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved  


Epsilon=10^(-5);
H=zeros(size(z,1),size(z,2));
idx1=find(z>Epsilon);
idx2=find(z<Epsilon & z>-Epsilon);
H(idx1)=1;
for i=1:length(idx2)
    H(idx2(i))=1/2*(1+z(idx2(i))/Epsilon+1/pi*sin(pi*z(idx2(i))/Epsilon));
end
end
function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);
end
function KG = kappa(I)
% get curvature information of input image
% input: 2D image I
% output: curvature matrix KG

% Copyright (c) 2009, 
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved  

I = double(I);
[m,n] = size(I);
P = padarray(I,[1,1],1,'pre');
P = padarray(P,[1,1],1,'post');

% central difference
fy = P(3:end,2:n+1)-P(1:m,2:n+1);
fx = P(2:m+1,3:end)-P(2:m+1,1:n);
fyy = P(3:end,2:n+1)+P(1:m,2:n+1)-2*I;
fxx = P(2:m+1,3:end)+P(2:m+1,1:n)-2*I;
fxy = 0.25.*(P(3:end,3:end)-P(1:m,3:end)+P(3:end,1:n)-P(1:m,1:n));
G = (fx.^2+fy.^2).^(0.5);
K = (fxx.*fy.^2-2*fxy.*fx.*fy+fyy.*fx.^2)./((fx.^2+fy.^2+eps).^(1.5));
KG = K.*G;
KG(1,:) = eps;
KG(end,:) = eps;
KG(:,1) = eps;
KG(:,end) = eps;
KG = KG./max(max(abs(KG)));
end

function indicator = checkstop(old,new,dt)
% indicate whether we should performance further iteraions or stop

% Copyright (c) 2009, 
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved  

    ind = find(abs(new)<=.5);
    M = length(ind);
    Q = sum(abs(new(ind)-old(ind)))./M;
    if Q<=dt*.18^2
        indicator = 1;
    else
        indicator = 0;
    end

return
end
