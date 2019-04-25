clear; close all; clc
% Video analysis on Loann's data


%% Loading video:

% Charge video:
vid = VideoReader('/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish6/Light On/Films/Fish On_6_8.avi');
vid = VideoReader('/home/ljp/Science/Loann/Free-Swimming-Postural-Stabilisation/no_feedback/Mesure/Fish9/5.avi');
disp(vid);
hv = vid.Height;
wv = vid.Width;
nv = vid.NumberofFrames;


%% Basic analysis:

% Animated video from images:
% figure; for i = 1:nv; temp = read(vid, i); imshow(temp); pause(1/vid.FrameRate); end

% % 3D-matrix with whole video:
% full_vid = zeros(hv, wv, nv);
% tic
% lprint = 0;
% for i = 1:nv
%     full_vid(:, :, i) = read(vid, i);
%     if mod(i, 100) == 0
%         fprintf(repmat('\b', 1, lprint));
%         lprint = fprintf('Iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, nv, toc]);
%     end
% end
% % Variance associated:
% var_vid = var(full_vid, [], 3);

% % Distances between pixels values:
% dist_vid = zeros(hv, wv);
% frame1 = read(vid, 1);
% tic; lprint = 0;
% figure
% for i = 2:nv
%     frame2 = read(vid, i);
%     frame_temp = double(abs(frame2 - frame1));
%     frame_temp = conv2(frame_temp, ones(3, 3)/9, 'same');
%     frame_temp(frame_temp < 5) = 0;
%     frame_temp(frame_temp >= 5) = 1;
%     dist_vid = dist_vid + frame_temp;
%     frame1 = frame2;
%     if mod(i, 100) == 0
%         fprintf(repmat('\b', 1, lprint));
%         lprint = fprintf('Iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, nv, toc]);
%         image(frame_temp, 'CDataMapping', 'scaled'); colorbar
%     end
% end
% dist_vid = dist_vid / (nv - 1);
% figure; image(dist_vid, 'CDataMapping', 'scaled'); colorbar


%% Correcting flashes of light:

% mean_vid = zeros(nv, 1);
% tic; lprint = 0;
% for i = 1:nv
%     frame_temp = read(vid, i);
%     mean_vid(i) = mean(frame_temp(:));
%     if mod(i, 100) == 0
%         fprintf(repmat('\b', 1, lprint));
%         lprint = fprintf('Iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, nv, toc]);
%     end
% end

% Inserting mean in previous algorithm:
% Parameters:
frame_lim = 7;
lconv = 3;
% Distances between pixels values:
mean_vid = zeros(nv, 1);
dist_vid = zeros(hv, wv);
frame1 = double(read(vid, 1));
mean_base = mean(frame1(:));
mean_vid(1) = mean_base;
tic; lprint = 0;
figure
for i = 2:nv
    % Load frame:
    frame2 = double(read(vid, i));
    % Subtract and save mean:
    frame2b = frame2 - mean(frame2(:)) + mean_base;
    mean_vid(i) = mean(frame2(:));
    % Compute difference with former frame:
    frame_temp = abs(frame2b - frame1);
    % Convolve and binarize:
    frame_temp = conv2(frame_temp, ones(lconv, lconv)/lconv^2, 'same');
    frame_temp(frame_temp < frame_lim) = 0;
    frame_temp(frame_temp >= frame_lim) = 1;
    % Add and change frame:
    dist_vid = dist_vid + frame_temp;
    frame1 = frame2b;
    % Indication:
    if mod(i, 50) == 0
        fprintf(repmat('\b', 1, lprint));
        lprint = fprintf('Iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, nv, toc]);
        imshow((frame_temp == 1) .* frame2)
    end
end
dist_vid = dist_vid / (nv - 1);
figure; image(dist_vid, 'CDataMapping', 'scaled'); colorbar


%% For Louis's movies:

% Parameters:
frame_lim = 5;
lconv = 1;
% Distances between pixels values:
mean_vid = zeros(nv, 1);
dist_vid = zeros(hv, wv);
frame1 = double(read(vid, 1));
frame1 = frame1(:, :, 1);
mean_base = mean(frame1(:));
mean_vid(1) = mean_base;
tic; lprint = 0;
figure
for i = 2:nv
    % Load frame:
    frame2 = double(read(vid, i));
    frame2 = frame2(:, :, 1);
    % Subtract and save mean:
    frame2b = frame2 - mean(frame2(:)) + mean_base;
    mean_vid(i) = mean(frame2(:));
    % Compute difference with former frame:
    frame_temp = abs(frame2b - frame1);
    % Convolve and binarize:
    frame_temp = conv2(frame_temp, ones(lconv, lconv)/lconv^2, 'same');
    frame_temp(frame_temp < frame_lim) = 0;
    frame_temp(frame_temp >= frame_lim) = 1;
    % Add and change frame:
    dist_vid = dist_vid + frame_temp;
    frame1 = frame2b;
    % Indication:
    if mod(i, 1) == 0
        fprintf(repmat('\b', 1, lprint));
        lprint = fprintf('Iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, nv, toc]);
        image((frame_temp == 1) .* frame2, 'CDataMapping', 'scaled'); colorbar
    end
end
dist_vid = dist_vid / (nv - 1);
figure; image(dist_vid, 'CDataMapping', 'scaled'); colorbar


%% For Louis's movies, using autocorrelation:

% Working with pooling so that it is easier:
pool_param = 20;
pool_mat = zeros(ceil(hv/pool_param), ceil(wv/pool_param), nv);
tic; lprint = 0;
for i = 1:nv
    % Load frame:
    frame = double(read(vid, i));
    frame = frame(:, :, 1);
    % Pool:
    frame = a02_avrg_pooling(frame, pool_param);
    % Insert in pool_mat:
    pool_mat(:, :, i) = frame;
    % Indication:
    if mod(i, 100) == 0
        fprintf(repmat('\b', 1, lprint));
        lprint = fprintf('For pooling, iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, nv, toc]);
    end
end

% Finding number of images per period:
corr_vid = zeros(nv, 1);
frame1 = pool_mat(:, :, 1);
mean_vid1 = mean(frame1(:));
tic; lprint = 0;
for i = 2:nv
    % Load frame:
    frame2 = pool_mat(:, :, i);
    % Subtract and save mean:
    frame2 = frame2 - mean(frame2(:)) + mean_vid1;
    % Add mean difference inverse:
    frame_temp = abs(frame2 - frame1);
    corr_vid(i) = 1 ./ mean(frame_temp(:));
    % Indication:
    if mod(i, 100) == 0
        fprintf(repmat('\b', 1, lprint));
        lprint = fprintf('For autocorrelation, iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, nv, toc]);
    end
end
% Finding rotation period:
[period_val, period_ind] = findpeaks(corr_vid);
temp_ind = period_ind(period_val > period_val(1)); % Finding high peaks
period = temp_ind(1) - 1;
% Now computing mean background for all points in the period:
background = zeros(hv, wv, period);
tic; lprint = 0;
for i = 1:period
    % Attributing right pool matrix:
    framei = pool_mat(:, :, i);
    framei = framei - mean(framei(:)) + mean_vid1;
    % Correlation vector:
    corr_vidi = zeros(nv, 1);
    for j = 1:nv
        % Load frame:
        framej = pool_mat(:, :, j);
        % Subtract and save mean:
        framej = framej - mean(framej(:)) + mean_vid1;
        % Add mean difference inverse:
        frame_temp = abs(framej - framei);
        corr_vidi(j) = 1 ./ mean(frame_temp(:));
    end
    corr_vidi(i) = 0;
    % Finding rotation period:
    [period_val, period_ind] = findpeaks(corr_vidi);
    if i == 1
        temp_ind = period_ind(period_val > corr_vidi(i+1)); % Finding high peaks
    elseif i == period
        temp_ind = period_ind(period_val > corr_vidi(i-1)); % Finding high peaks
    else
        temp_ind = period_ind(period_val > min([corr_vidi(i-1), corr_vidi(i+1)])); % Finding high peaks
    end
    temp_ind = [i; temp_ind]; % Adding initial value
    ltind = length(temp_ind);
    % Adding each frame one after the other to background:
    for k = 1:ltind
        % Loading frame:
        frame_temp = double(read(vid, temp_ind(k)));
        frame_temp = frame_temp(:, :, 1);
        frame_temp = frame_temp - mean(frame_temp(:)) + mean_vid1;
        % Add to background:
        background(:, :, i) = background(:, :, i) + frame_temp/ltind;
    end 
    % Indication:
    if mod(i, 10) == 0
        fprintf(repmat('\b', 1, lprint));
        lprint = fprintf('For autocorrelation, iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, period, toc]);
    end
end

figure
for i = 1:nv
    % Loading frame:
    frame_temp = double(read(vid, i));
    frame_temp = frame_temp(:, :, 1);
    frame_temp = frame_temp - mean(frame_temp(:)) + mean_vid1;
    % Compute absolute difference with background:
    backind = mod(i-1, period) + 1;
    diff_temp = abs(frame_temp - background(:, :, backind));
    % Plotting absolute difference:
    pause(0.001)
    image(diff_temp, 'CDataMapping', 'scaled')
end
    


%% New idea, autocorrelation:
% Using correlation to first frame, we are going to be able to tell after
% how long the tank has done a whole turn, taking maximum correlation. Then
% we will compare pixels from first frame to period frame, and find
% redundant pixels on the side.

% Distances between pixels values:
corr_vid = zeros(nv, 1);
frame1 = double(read(vid, 201));
tic; lprint = 0;
for i = 2:nv
    % Load frame:
    frame2 = double(read(vid, i));
    % Subtract and save mean:
    frame2 = frame2 - mean_vid(i) + mean_vid(1);
    % Add mean difference:
    frame_temp = abs(frame2 - frame1);
    corr_vid(i) = 1 ./ mean(frame_temp(:));
    % Indication:
    if mod(i, 100) == 0
        fprintf(repmat('\b', 1, lprint));
        lprint = fprintf('Iteration %.0f out of %.0f, done in %.3f seconds. \n', [i, nv, toc]);
    end
end
% Alright, period is approx. 1330 frames.


%% We are going to build a database, comprising points for:
%  - End of tail
%  - Middle of tail
%  - Beginning of tail
%  - Each eye
%  - Beginning of the fish (that is the nose)
%  For a total of 6 points per frame, respectively in previous order.

% Trial to see if easy to get points:
numex = 3;
pdata = zeros(2, 6, numex);

for i = 1:numex
    % Plot image and allow zoom:
    h.myfig = figure('units','normalized','outerposition',[0 0 1 1]);
    frame_temp = read(vid, i);
    if i ~= 1
        frame_temp(pdata(1, :, i-1), pdata(2, :, i-1)) = 0;
    end
    imshow(frame_temp)
    zoom on
    while waitforbuttonpress 
    end
    zoom off
    % Loop over the 6 points needed:
    for ipt = 1:6
        cursorobj = datacursormode(h.myfig);
        cursorobj.SnapToDataVertex = 'on'; % Snap to our plotted data, on by default
        cursorobj.Enable = 'on'; % Turn on the data cursor, hold alt to select multiple points
        while ~waitforbuttonpress 
            % waitforbuttonpress returns 0 with click, 1 with key press
            % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
        end
        cursorobj.Enable = 'off';
        % Add information in pdata:
        mypoints = getCursorInfo(cursorobj);
        pdata(1, ipt, i) = mypoints.Position(1);
        pdata(2, ipt, i) = mypoints.Position(2);
    end
    close(h.myfig)
end



tg = a02_avrg_pooling(temp1, 10);








