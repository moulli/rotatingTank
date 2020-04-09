clear; close all; clc
addpath(genpath('/home/ljp/Science/balancoire_data/balancoire_matlab_files'))
addpath(genpath('/home/ljp/Science/Hippolyte/rotatingTank'))


%% Prepare output of file

washingOut = struct;


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
numframes = 255; %floor(vid.Duration * vid.FrameRate);
width = vid.Width;
height = vid.Height;

% Add to output structure
washingOut.path = vidpath;


%% Get background and masks

[background, masks] = getMasks(vidpath, 'numframes', numframes);

% Add to output structure
washingOut.preprocessing.background = background;
washingOut.preprocessing.masks = masks;


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
        outmask = 8; % arbitrary number for first mask
    else
        outmask = choseMask(im, background, masks, dp(frame-1, :), max_diff_disp);
    end
    maskused(frame) = outmask;
    
    % Get highest pixel
    new_im = (background - im) .* masks(:, :, outmask);
    [dpframex, dpframey] = find(new_im == max(new_im(:)));
    newdp = [dpframex(1), dpframey(1)];
    dp(frame, :) = newdp;
    % Get wider mask, in case tail is outside ROI
    outmask = min([outmask+2, size(masks, 3)]);
    new_im = (background - im) .* masks(:, :, outmask);
    
    % Make a mask to isolate neighbourhood
    [xnei, ynei] = meshgrid(-(newdp(2)-1):(width-newdp(2)), -(newdp(1)-1):(height-newdp(1))); 
    masknei = ((xnei.^2 + ynei.^2) <= rneighb^2);
    new_im = masknei .* new_im;
    
    % Pad new_im
    new_im = padarray(new_im, [rneighb, rneighb], 0, 'both');
    % Save neighbourhood
    fish_neighb(:, :, frame) = new_im(newdp(1):newdp(1)+2*rneighb, newdp(2):newdp(2)+2*rneighb);
    
end

% POST TREATMENT: x axis is y axis on matlab
dp = fliplr(dp);

% Add to output structure
washingOut.images.parameters.max_diff_disp = max_diff_disp;
washingOut.images.parameters.rneighb = rneighb;
washingOut.images.dp = dp;
washingOut.images.maskused = maskused;
washingOut.images.fish_neighb = fish_neighb;


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
    
end

% POST TREATMENT: x axis is y axis on matlab
for frame = 1:numframes
    fpts{frame} = fliplr(fpts{frame});
end

% POST TREATMENT: remove anomaly points
% If orientation for x-coordinate and y-coordinate is inverted for only one
% time point, we will say this time point is the mean of the former and the
% later position
modified = zeros(numframes, 1);
for frame = 2:(numframes-1)
    % Orientation from before
    ori_m1 = [fpts{frame-1}(2, 1)-fpts{frame-1}(1, 1), fpts{frame-1}(2, 2)-fpts{frame-1}(1, 2)];
    ori_m1 = (ori_m1 > 0);
    % Orientation from now
    ori_n = [fpts{frame}(2, 1)-fpts{frame}(1, 1), fpts{frame}(2, 2)-fpts{frame}(1, 2)];
    ori_n = (ori_n < 0);
    % Orientation from before
    ori_p1 = [fpts{frame+1}(2, 1)-fpts{frame+1}(1, 1), fpts{frame+1}(2, 2)-fpts{frame+1}(1, 2)];
    ori_p1 = (ori_p1 > 0);
    % Correct if necessary
    if isequal(ori_m1, ori_n) && isequal(ori_n, ori_p1)
        modified(frame) = 1;
        fpts{frame}(2:end, :) = (fpts{frame-1}(2:end, :) + fpts{frame+1}(2:end, :)) / 2;
    end    
end

% Get info from fpts
ptss = zeros(num_segments+1, 2, numframes);
for frame = 1:numframes
    for pt = 1:num_segments+1
        for coord = 1:2
            ptss(pt, coord, frame) = fpts{frame}(pt, coord) + dp(frame, coord);
            % NB: equivalent to ptss(:, :, frame) = fpts{frame} + dp(frame, :);
        end
    end
end

% Add to output structure
washingOut.tracking.parameters.num_segments = num_segments;
washingOut.tracking.parameters.tail_length = tail_length;
washingOut.tracking.parameters.body_length = body_length;
washingOut.tracking.parameters.inertia = inertia;
washingOut.tracking.parameters.num_pix1 = num_pix1;
washingOut.tracking.parameters.num_pix2 = num_pix2;
washingOut.tracking.parameters.initial_box = initial_box;
washingOut.tracking.fpts = fpts;
washingOut.tracking.modified = modified;
washingOut.tracking.ptss = ptss;


%% Plot fish position

% Animation
figure
inst = 1;
for frame = inst:numframes
    % Plot
    hold on
    image(washingOut.images.fish_neighb(:, :, frame), 'CDataMapping', 'scaled')
    plot(washingOut.tracking.fpts{frame}(:, 1), washingOut.tracking.fpts{frame}(:, 2), 'r')
    axis equal
    pause(0.05)
    hold off
end


%% Plot displacements

figure
dt = 1 / vid.FrameRate;
time = (0:1:(numframes-1)) * dt;
subplot(2, 1, 1)
hold on
for i = 1:num_segments+1
    tseriex = reshape(washingOut.tracking.ptss(i, 1, :), [numframes, 1]);
    plot(time, tseriex)
end
plot(time, washingOut.tracking.modified*200)
title('x-axis displacement', 'Interpreter', 'latex')
xlabel('time [s]', 'Interpreter', 'latex')
ylabel('x-coordinate', 'Interpreter', 'latex')
subplot(2, 1, 2)
hold on
for i = 1:num_segments+1
    tseriey = reshape(washingOut.tracking.ptss(i, 2, :), [numframes, 1]);
    plot(time, tseriey)
end
plot(time, washingOut.tracking.modified*200)
title('y-axis displacement', 'Interpreter', 'latex')
xlabel('time [s]', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')


%% Identifying fish rotation along rostro-caudal axis

% Negative log-likelihood for gaussian fitting
negll = zeros(numframes, 2);
mu2 = zeros(2, 2, numframes);
mu3 = zeros(3, 2, numframes);

% Loop over all frames
for frame = 1:numframes
    
    % Read frame
    new_im = fish_neighb(:, :, frame);
    
    % Get highest pixels
    num_pts = 91*91;
    f_ind1 = find(new_im >= quantile(new_im(:), 1 - num_pts/numel(new_im)));
    pix_value1 = new_im(f_ind1);
    [fx1, fy1] = ind2sub(size(new_im), f_ind1);
    
    % Make dataset for gaussian fitting
    gauss_set = zeros(0, 2);
    for gpt = 1:length(pix_value1)
        gauss_set = cat(1, gauss_set, [fx1(gpt), fy1(gpt)].*ones(round(pix_value1(gpt)), 2));
    end
    gauss_set = gauss_set + 0.1*(rand(size(gauss_set)) - 0.5); % add noise (not mandatory)
    
    %
    gm2 = fitgmdist(gauss_set, 2);
    gm3 = fitgmdist(gauss_set, 3);
    negll(frame, :) = [gm2.NegativeLogLikelihood, gm3.NegativeLogLikelihood];
    mu2(:, :, frame) = gm2.mu;
    mu3(:, :, frame) = gm3.mu;
%     mu = gm.mu;
%     disteye = sqrt((mu(:, 1) - mu(:, 1)').^2 + (mu(:, 2) - mu(:, 2)').^2);
%     inter_eye(frame) = mean(mean(disteye));
   
end

% Comparing area between 3 points
mu2 = sort(mu2);
mu3 = sort(mu3);
tarea = zeros(numframes, 1);
for frame = 1:numframes
    temp = mu3([1, 2, 3, 1], :, frame);
    tarea(frame) = polyarea(temp(:, 1), temp(:, 2));
end


    
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

% Comparing distribution of pixels
distf = zeros(255, numframes);
for frame = 1:numframes
    for intens = 1:255
        fnei = fish_neighb(:, :, frame);
        distf(intens, frame) = sum(intens-1 <= fnei(:) & fnei(:) < intens);
    end
end
distf = distf ./ sum(distf);
figure
hold on
for frame = 1:numframes
    plot(3:255, distf(3:end, frame))
end
abc = zeros(3, numframes);
for frame = 1:numframes
    gm = fit((3:255)', distf(3:end, frame), 'gauss1');
    abc(:, frame) = [gm.a1; gm.b1; gm.c1];
end
figure
hold on
for i = 1:3
    plot(abc(i, :))
end

% Maximul pooling
spool = 3;
sim = 2*rneighb+1;
outi = zeros(floor(sim/spool));
for i = 1:floor(sim/spool)
    for j = 1:floor(sim/spool)
        outemp = new_im(1+spool*(i-1):spool*i, 1+spool*(j-1):spool*j);
        outi(i, j) = max(outemp(:));
    end
end
figure
image(outi, 'CDataMapping', 'scaled')
axis equal


    





