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
numframes = 200; %floor(vid.Duration * vid.FrameRate);
width = vid.Width;
height = vid.Height;


%% Get background and masks

[background, masks] = getMasks(vidpath, 'numframes', numframes);


%% Loop over video to get fish position

% Redefine video reader
vid = VideoReader(vidpath);

% Define maximum displacement differential to switch mask in needed
max_diff_disp = 20;
% Define neighbourhood to consider around fish
rneighb = 45;
% Segment tracking parameters
num_segments = 5;
tail_length = 30;
body_length = 0;
inertia = 0.5;
num_pix1 = 10;
num_pix2 = 150;
initial_box = 0.6;

% Define outputs
dp = zeros(numframes, 2);
fpts = cell(numframes, 1);
maskused = zeros(numframes, 1);

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
    
    % Identifying important points
    new_im = (255 - new_im);
    segment_pts = segmentTracking_WM(new_im, 'num_segments', num_segments, ...
                        'tail_length', tail_length, 'body_length', body_length, ...
                        'inertia', inertia, 'num_pix1', num_pix1, ...
                        'num_pix2', num_pix2, 'initial_box', initial_box);
    fpts{frame} = segment_pts;
    
%     % Identify fish rotation along rostro-caudal axis
%     num_pts = 10;
%     f_ind1 = find(new_im <= quantile(new_im(:), num_pts/numel(new_im)));
%     pix_value1 = (256 - double(new_im(f_ind1)));
%     com1(frame, :) = pix_value1' * [fx1, fy1] ./ sum(pix_value1);
%     inter_eye(frame) = sqrt(sum((com1(frame, :) - newdp).^2));
    
%     % Identify fish rotation along rostro-caudal axis
%     num_pts = 10;
%     f_ind1 = find(new_im <= quantile(new_im(:), num_pts/numel(new_im)));
%     [fx1, fy1] = ind2sub(size(new_im), f_ind1);
%     dist_eye = sqrt((fx1 - fx1').^2 + (fy1 - fy1').^2);
%     inter_eye(frame) = mean(mean(dist_eye));

    % Identify fish rotation along rostro-caudal axis
    num_pts = 50;
    f_ind1 = find(new_im <= quantile(new_im(:), num_pts/numel(new_im)));
    pix_value1 = (256 - (new_im(f_ind1)));
    [fx1, fy1] = ind2sub(size(new_im), f_ind1);
    gauss_set = zeros(0, 2);
    for gpt = 1:length(pix_value1)
        gauss_set = cat(1, gauss_set, [fx1(gpt), fy1(gpt)].*ones(round(pix_value1(gpt)), 2));
    end
    gauss_set = gauss_set + 0.01*(rand(size(gauss_set)) - 0.5);
    gm = fitgmdist(gauss_set, 3);
    inter_eye(frame) = gm.NegativeLogLikelihood;
%     mu = gm.mu;
%     disteye = sqrt((mu(:, 1) - mu(:, 1)').^2 + (mu(:, 2) - mu(:, 2)').^2);
%     inter_eye(frame) = mean(mean(disteye));
end

% Post treatment (x axis is y axis on matlab)
dp = fliplr(dp);
for frame = 1:numframes
    fpts{frame} = fliplr(fpts{frame});
end


%% Plot fish position

% Redefine video reader
vid = VideoReader(vidpath);

% Get info from fpts
ptss = zeros(2*(num_segments+1), numframes);
for i = 1:numframes
    for j = 1:num_segments+1
        for k = 1:2
            ptss(k+2*(j-1), i) = fpts{i}(j, k);
        end
    end
end

% Animation
figure
inst = 1;
for i = 1:(inst-1)
    im = readFrame(vid);
end
for frame = inst:numframes
    im = readFrame(vid);
    im = flipud(double(im));
    % Treat image
    new_im = (background - im) .* masks(:, :, maskused(frame));
    [xnei, ynei] = meshgrid(-(dp(frame, 1)-1):(width-dp(frame, 1)), -(dp(frame, 2)-1):(height-dp(frame, 2))); 
    masknei = ((xnei.^2 + ynei.^2) <= rneighb^2);
    new_im = 255 - masknei .* new_im;
    % Plot
    hold on
    image(new_im, 'CDataMapping', 'scaled')
    plot(ptss(1:2:end, frame), ptss(2:2:end, frame), 'r')
    axis equal
    pause(0.01)
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

figure
subplot(2, 1, 1)
hold on
for i = 1:num_segments+1
    plot(diff(ptss(1 + 2*(i-1), :)))
end
subplot(2, 1, 2)
hold on
for i = 1:num_segments+1
    plot(diff(ptss(2*i, :)))
end

figure
hold on
plot(diff(ptss(end-1, :)-ptss(1, :)))
plot(diff(ptss(end, :)-ptss(2, :)))

figure
ptss1 = sqrt(ptss(1, :).^2 + ptss(2, :).^2);
ptss2 = sqrt(ptss(end-1, :).^2 + ptss(end, :).^2);
plot(diff(ptss2 - ptss1))


%% Post treatment of data
% So I don't think there is a simple possibility of having all data good
% when tracking, therefore I am going to do some post treatment and replace
% wrong tracking data with nan values

angles = zeros(num_segments, numframes);
for i = 1:num_segments
    for j = 1:numframes
        xi = ptss(1 + 2*i, j) - ptss(1, j);
        yi = ptss(2*(i+1), j) - ptss(2, j);
        angles(i, j) = atan(yi / xi);
    end
end
figure
hold on
for i = 1:num_segments
    plot((angles(i, :)))    
end

figure
subplot(3, 1, 1)
hold on
for i = 2:num_segments+1
    plot(ptss(1 + 2*(i-1), :) - ptss(1, :))
end
subplot(3, 1, 2)
hold on
for i = 2:num_segments+1
    plot(ptss(2*i, :) - ptss(2, :))
end
subplot(3, 1, 3)
hold on
for i = 2:num_segments+1
    plot(sqrt(diff(ptss(1 + 2*(i-1), :) - ptss(1, :)).^2 + diff(ptss(2*i, :) - ptss(2, :)).^2))
end

% After reflexion, let's keep these outliers. They represent only a
% minority of points, and usually happen during a bout, so it is not the
% worst scenario...





