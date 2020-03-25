clear; close all; clc
addpath(genpath('/home/ljp/Science/balancoire_data/balancoire_matlab_files'))
addpath(genpath('/home/ljp/Science/Hippolyte/rotatingTank'))


%% Get video

vidnum = 3;
if vidnum == 1
    vidpath = '/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish1/Light On/Films/Fish On_1_12.avi';
elseif vidnum == 2
    vidpath = '/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish13/Light Off/Films/Fish Off_13_6.avi';
elseif vidnum == 3
    vidpath = '/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish18/Light On/Films/Fish On_18_0.avi';
end
vid = VideoReader(vidpath);
numframes = floor(vid.Duration * vid.FrameRate);
width = vid.Width;
height = vid.Height;


%% Get background and masks

[background, masks] = getMasks(vidpath, 'numframes', numframes);


%% Loop over video and save neighbourhood of fish

% Redefine video reader
vid = VideoReader(vidpath);

% Define maximum displacement differential to switch mask in needed
max_diff_disp = 20;
% Define neighbourhood to consider around fish
rneighb = 45;

% Define outputs
dp = zeros(numframes, 2);
maskused = zeros(numframes, 1);

% Fish neighbourhood
fish_neighb = zeros(2*rneighb+1, 2*rneighb+1, numframes);

% Loop over all frames
for frame = 1:numframes
    
    % Read frame
    im = readFrame(vid);
    im = flipud(double(im));
    
    % Chose right mask
    if frame == 1
        outmask = 10; % arbitrary number for first mask
    else
        outmask = choseMask(im, background, masks, dp(frame-1, :), max_diff_disp);
    end
    maskused(frame) = outmask;
    
    % Get highest pixel
    new_im = (background - im) .* masks(:, :, outmask);
    [dpframex, dpframey] = find(new_im == max(new_im(:)));
    newdp = [dpframex(1), dpframey(1)];
    dp(frame, :) = newdp;
    
    % Make a mask to isolate neighbourhood
    [xnei, ynei] = meshgrid(-(newdp(2)-1):(width-newdp(2)), -(newdp(1)-1):(height-newdp(1))); 
    masknei = ((xnei.^2 + ynei.^2) <= rneighb^2);
    new_im = masknei .* new_im;
    
    % Pad new_im
    new_im = padarray(new_im, [rneighb, rneighb], 0, 'both');
    % Save neighbourhood
    fish_neighb(:, :, frame) = new_im(newdp(1):newdp(1)+2*rneighb, newdp(2):newdp(2)+2*rneighb);
    
end

% Post treatment (x axis is y axis on matlab)
dp = fliplr(dp);


%% Loop over video and save neighbourhood of fish

% Segment tracking parameters
num_segments = 5;
tail_length = 30;
body_length = 0;
inertia = 0.5;
num_pix1 = 10;
num_pix2 = 150;
initial_box = 0.6;

% Define outputs
fpts = cell(numframes, 1);

% Loop over all frames
for frame = 1:numframes
    
    % Read frame
    new_im = fish_neighb(:, :, frame);
    
    % Identifying important points
    new_im = (255 - new_im);
    segment_pts = segmentTracking_WM(new_im, 'num_segments', num_segments, ...
                        'tail_length', tail_length, 'body_length', body_length, ...
                        'inertia', inertia, 'num_pix1', num_pix1, ...
                        'num_pix2', num_pix2, 'initial_box', initial_box);
    fpts{frame} = segment_pts;
    
% %     % Identify fish rotation along rostro-caudal axis
% %     num_pts = 10;
% %     f_ind1 = find(new_im <= quantile(new_im(:), num_pts/numel(new_im)));
% %     pix_value1 = (256 - double(new_im(f_ind1)));
% %     com1(frame, :) = pix_value1' * [fx1, fy1] ./ sum(pix_value1);
% %     inter_eye(frame) = sqrt(sum((com1(frame, :) - newdp).^2));
%     
% %     % Identify fish rotation along rostro-caudal axis
% %     num_pts = 10;
% %     f_ind1 = find(new_im <= quantile(new_im(:), num_pts/numel(new_im)));
% %     [fx1, fy1] = ind2sub(size(new_im), f_ind1);
% %     dist_eye = sqrt((fx1 - fx1').^2 + (fy1 - fy1').^2);
% %     inter_eye(frame) = mean(mean(dist_eye));
% 
%     % Identify fish rotation along rostro-caudal axis
%     num_pts = 50;
%     f_ind1 = find(new_im <= quantile(new_im(:), num_pts/numel(new_im)));
%     pix_value1 = (256 - (new_im(f_ind1)));
%     [fx1, fy1] = ind2sub(size(new_im), f_ind1);
%     gauss_set = zeros(0, 2);
%     for gpt = 1:length(pix_value1)
%         gauss_set = cat(1, gauss_set, [fx1(gpt), fy1(gpt)].*ones(round(pix_value1(gpt)), 2));
%     end
%     gauss_set = gauss_set + 0.01*(rand(size(gauss_set)) - 0.5);
%     gm = fitgmdist(gauss_set, 3);
%     inter_eye(frame) = gm.NegativeLogLikelihood;
% %     mu = gm.mu;
% %     disteye = sqrt((mu(:, 1) - mu(:, 1)').^2 + (mu(:, 2) - mu(:, 2)').^2);
% %     inter_eye(frame) = mean(mean(disteye));
end

% Post treatment (x axis is y axis on matlab)
for frame = 1:numframes
    fpts{frame} = fliplr(fpts{frame});
end

% Get info from fpts
ptss = zeros(2*(num_segments+1), numframes);
for i = 1:numframes
    for j = 1:num_segments+1
        for k = 1:2
            ptss(k+2*(j-1), i) = fpts{i}(j, k) + dp(i, k);
        end
    end
end


%% Plot fish position

% Redefine video reader
vid = VideoReader(vidpath);

% Animation
figure
inst = 1;
for frame = inst:numframes
    % Plot
    hold on
    image(fish_neighb(:, :, frame), 'CDataMapping', 'scaled')
    plot(fpts{frame}(:, 1), fpts{frame}(:, 2), 'r')
    axis equal
    pause(0.05)
    hold off
end


%% Plot displacements

figure
subplot(2, 1, 1)
hold on
for i = 1:num_segments+1
    plot(ptss(1 + 2*(i-1), :))
end
subplot(2, 1, 2)
hold on
for i = 1:num_segments+1
    plot(ptss(2*i, :))
end





