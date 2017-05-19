function redraw(frame, vidObj,mu,a_1,phi,timestep,epsilon,K,num_inter,K_phi,c0)
% REDRAW  Process a particular frame of the video
%   REDRAW(FRAME, VIDOBJ)
%       frame  - frame number to process
%       vidObj - VideoReader object

% Read frame
f = vidObj.read(frame);

% Get edge
f2 = edge(rgb2gray(f), 'canny');

% Overlay edge on original image
f3 = bsxfun(@plus, f,  uint8(255*f2));

% Display
image(f3); axis image off

end
