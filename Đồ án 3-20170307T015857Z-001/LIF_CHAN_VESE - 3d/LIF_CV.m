function [phi] = LIV_CV(sigma_phi,num_inter,sigma,file_name,mu,a_1,timestep,epsilon)
 
c0 =2;
Img = load(file_name);
Img = Img.B;
%sigma =3;% the key parameter which needs to be tuned properly.
%sigma_phi = 0.45;% the variance of regularized Gaussian kernel
K = fspecial3('gaussian',round(sigma_phi*(4*sqrt(2*log(2)))));
K_phi = fspecial3('gaussian',round(sigma_phi*(4*sqrt(2*log(2)))));
I=Img(:,:,1);
[nrow,ncol,nhigh] = size(Img);
margin = 1; 
phi = zeros(size(Img)); 
phi(margin:end-margin, margin:end-margin, margin:end-margin) = 1; 
phi = (phi-.5); 
   %Draw mask
      % figure;imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;axis equal;
      % text(6,6,'Left click to get points, right click to get end point','FontSize',[12],'Color', 'g');
      % BW=roipoly;
       %phi=c0*2*(0.5-BW);
%        phi=zeros(nrow,ncol,nhigh);
%        phi=phi-1;
%        phi(nrow,:,:)=0;
%        phi(:,ncol,:)=0;
%        phi(:,:,nhigh)=0;
       %hold on;
      % [c,h]=contour(phi,[0 0],'r');
      %set(h, 'linewidth', 2.5);
      % hold off;
       %end draw map;

%pause(0.01);
%figure;  imagesc(Img,[0 255]);colormap(gray);hold on;
%contour(phi,[0 0],'b');

%timestep = 0.005;
%epsilon = 1.0;

for n = 1:num_inter   
   old=phi;
      [force,phi,f1,f2,Hphi]= Caculate_next_phi(mu,a_1,Img,phi,timestep,epsilon,sigma);
         
       phi = imgaussfilt3(phi,sigma_phi);
       new=phi; 
       iso = isosurface(phi);
       if exist('h','var') && all(ishandle(h)), delete(h); end
    h = patch(iso,'facecolor','w');  axis equal;  view(3); 
    set(gcf,'name', sprintf('#iters = %d',n));
    drawnow; 
     % indicator=checkstop(old,new,timestep);
      %if(indicator) 
         % imagesc(Img,[0 255]);colormap(gray)
     % hold on;[c,h]=contour(phi,[0 0],'b');
     % set(h, 'linewidth', 2.5);
      %iterNum=[num2str(n), ' iterations'];        
     % title(iterNum);        
   %   hold off;
     % figure;
     % mesh(phi);
     % return;
      %end
      if mod(n,40)==0
     % pause(0.0001);
      %imagesc(Img,[0 255]);colormap(gray)
     % hold on;[c,h]=contour(phi,[0 0],'b');
     % set(h, 'linewidth', 2.5);
      %iterNum=[num2str(n), ' iterations'];        
     % title(iterNum);        
     % hold off;
      end
end
figure;
slice = [10,15,20,25,30,35,40,45];
for j = 1:8
    subplot(2,4,j); imshow(Img(:,:,slice(j)),[]); hold on; 
    c = contours(phi(:,:,slice(j)),[0,0]);
    zy_plot_contours(c,'linewidth',2);
end
end
%figure;
%mesh(phi);
function [force,phi,f1,f2,Hphi] = Caculate_next_phi(mu,a_1,I,phi,timestep,epsilon,K)
inidx = find(phi>=0); % frontground index
outidx = find(phi<0); % background index
phi = NeumannBoundCond(phi);
P = double(I);
L = im2double(P(:,:,:)); % get one image component
force_image = 0; %
c1 = sum(sum(sum(L.*Heaviside1(phi))))/(length(inidx)+eps); % average inside of Phi0
c2 = sum(sum(sum(L.*(1-Heaviside1(phi)))))/(length(outidx)+eps); % verage outside of Phi0

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

f1 = imgaussfilt3(I.*H,K);
c1 = imgaussfilt3(H,K);
f1 = f1./c1;
f2 = imgaussfilt3(I.*(1-H),K);
c2 = imgaussfilt3(1-H,K);
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
H=zeros(size(z,1),size(z,2),size(z,3));
idx1=find(z>Epsilon);
idx2=find(z<Epsilon & z>-Epsilon);
H(idx1)=1;
for i=1:length(idx2)
    H(idx2(i))=1/2*(1+z(idx2(i))/Epsilon+1/pi*sin(pi*z(idx2(i))/Epsilon));
end;
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
function h = fspecial3(type,siz)
%FSPECIAL3 Create predefined 3-D filters.
%   H = FSPECIAL3(TYPE,SIZE) creates a 3-dimensional filter H of the
%   specified type and size. Possible values for TYPE are:
%
%     'average'    averaging filter
%     'ellipsoid'  ellipsoidal averaging filter
%     'gaussian'   Gaussian lowpass filter
%     'laplacian'  Laplacian operator
%     'log'        Laplacian of Gaussian filter
%
%   The default SIZE is [5 5 5]. If SIZE is a scalar then H is a 3D cubic
%   filter of dimension SIZE^3.
%
%   H = FSPECIAL3('average',SIZE) returns an averaging filter H of size
%   SIZE. SIZE can be a 3-element vector specifying the dimensions in
%   H or a scalar, in which case H is a cubic array.
%
%   H = FSPECIAL3('ellipsoid',SIZE) returns an ellipsoidal averaging filter.
%
%   H = FSPECIAL3('gaussian',SIZE) returns a centered Gaussian lowpass
%   filter of size SIZE with standard deviations defined as
%   SIZE/(4*sqrt(2*log(2))) so that FWHM equals half filter size
%   (http://en.wikipedia.org/wiki/FWHM). Such a FWHM-dependent standard
%   deviation yields a congruous Gaussian shape (what should be expected
%   for a Gaussian filter!).
%
%   H = FSPECIAL3('laplacian') returns a 3-by-3-by-3 filter approximating
%   the shape of the three-dimensional Laplacian operator. REMARK: the
%   shape of the Laplacian cannot be adjusted. An infinite number of 3D
%   Laplacian could be defined. If you know any simple formulation allowing
%   one to control the shape, please contact me.
%
%   H = FSPECIAL3('log',SIZE) returns a rotationally symmetric Laplacian of
%   Gaussian filter of size SIZE with standard deviation defined as
%   SIZE/(4*sqrt(2*log(2))).
%
%   Class Support
%   -------------
%   H is of class double.
%
%   Example
%   -------
%      I = single(rand(80,40,20));
%      h = fspecial3('gaussian',[9 5 3]); 
%      Inew = imfilter(I,h,'replicate');
%       
%   See also IMFILTER, FSPECIAL.
%
%   -- Damien Garcia -- 2007/08

narginchk(1,2)
type = lower(type);

if nargin==1
        siz = 5;
end

if numel(siz)==1
    siz = round(repmat(siz,1,3));
elseif numel(siz)~=3
    error('Number of elements in SIZ must be 1 or 3')
else
    siz = round(siz(:)');
end

switch type
    
    case 'average'
        h = ones(siz)/prod(siz);
        
    case 'gaussian'
        sig = siz/(4*sqrt(2*log(2)));
        siz   = (siz-1)/2;
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
        h = h/sum(h(:));
        
    case 'ellipsoid'
        R = siz/2;
        R(R==0) = 1;
        h = ones(siz);
        siz = (siz-1)/2;
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        I = (x.*x/R(1)^2+y.*y/R(2)^2+z.*z/R(3)^2)>1;
        h(I) = 0;
        h = h/sum(h(:));
        
    case 'laplacian'
        h = zeros(3,3,3);
        h(:,:,1) = [0 3 0;3 10 3;0 3 0];
        h(:,:,3) = h(:,:,1);
        h(:,:,2) = [3 10 3;10 -96 10;3 10 3];
        
    case 'log'
        sig = siz/(4*sqrt(2*log(2)));
        siz   = (siz-1)/2;
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
        h = h/sum(h(:));
        arg = (x.*x/sig(1)^4 + y.*y/sig(2)^4 + z.*z/sig(3)^4 - ...
            (1/sig(1)^2 + 1/sig(2)^2 + 1/sig(3)^2));
        h = arg.*h;
        h = h-sum(h(:))/prod(2*siz+1);
        
    otherwise
        error('Unknown filter type.')
        
end
end
        function u = ac_reinit(u)
% FUNCTION U_SDF = AC_REINIT(U)
% Reinitialization of the levelset U. Calculate the sign distance function
% of the zero-set of the input function U.
%   Input: 
%       U --- Input levelset. 
%   Output:
%       U_SDF --- Signed distance function based on the zero-set of the
%           input levelset.
%
% Example:
%     u = (imread('coins.png') > 100)-.5; 
%     u_sdf = ac_reinit(u); 
%     imshow(u_sdf,[]); colormap(jet); 


if ndims(u) == 2
    %% Extract zero-set.
    c = contours(u,[0,0]);
    xy = zy_extract_pt_from_contours(c);
    if isempty(xy), u = []; return; end
    %% Generate sign distance function.
    u0 = zeros(size(u));
    u0(sub2ind(size(u), round(xy(2,:)),round(xy(1,:)))) = 1;
    u = double(bwdist(u0)).*sign(u);
elseif ndims(u) == 3
%     %% Extract zero-set
%     iso = isosurface(u,0);
%     %% Generate sign distance function.
%     u0 = zeros(size(u));
%     idx = sub2ind(size(u), round(iso.vertices(:,1)), ...
%         round(iso.vertices(:,2)), round(iso.vertices(:,3)));
%     u0(idx) = 1;

    u0 = zy_binary_boundary_detection(uint8(u>0)); 
    % NOTE: bwdist for 3D is just too slow, so fast marching is applied.
%     u  = bwdist(u0).*sign(u); 
    u = ac_distance_transform_3d(u0).*sign(u); 
end
        end
        