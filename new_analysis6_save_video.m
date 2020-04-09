clear; close all; clc
folderpath = '/home/ljp/Science/rotatingTank_vidMatrices';

% fish parameters
f = 1;
l = ' On';
r = 12;

% load video
vidpath = strcat('/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish', ...
                              num2str(f), '/Light', l, '/Films/Fish' , l, '_', ...
                              num2str(f), '_', num2str(r), '.avi');
vid = VideoReader(vidpath);                          

% load structure
washpath = fullfile(folderpath, strcat('washingOut_', num2str(f), '_', l(2:end), '_', num2str(r), '.mat'));
load(washpath, 'washingOut');

% load image
im = readFrame(vid);
% im = flipud(im);

% get interesting data from washingOut
dp = washingOut.images.dp(1, :);
fpts = washingOut.tracking.fpts{1};
fpts_ = fpts - fpts(1, :);
tr = dp + fpts_;

% adapt data to image
dp = [dp(1), size(im, 1)-dp(2)];
tr = [tr(:, 1), size(im, 1)-tr(:, 2)];

% plot with tracking
figure
imshow(im)
hold on
scatter(dp(1), dp(2), 'y')
plot(tr(:, 1), tr(:, 2), 'y')


% video 
vid = VideoReader(vidpath);
numframes = floor(vid.Duration * vid.FrameRate);
% make video
vidname = fullfile(folderpath, strcat('vid_', num2str(f), '_', l(2:end), '_', num2str(r), '.avi'));
writerObj = VideoWriter(vidname);
writerObj.FrameRate = vid.FrameRate;
open(writerObj);
figure
for i = 1:numframes
    % Read frame and analyse it
    im = readFrame(vid);
    % get tracking
    dp = washingOut.images.dp(i, :);
    fpts = washingOut.tracking.fpts{i};
    fpts_ = fpts - fpts(1, :);
    tr = dp + fpts_;
    % adapt data to image
    dp = [dp(1), size(im, 1)-dp(2)];
    tr = [tr(:, 1), size(im, 1)-tr(:, 2)];
    % plot with tracking
    imshow(im)
    hold on
    scatter(dp(1), dp(2), 'y')
    plot(tr(:, 1), tr(:, 2), 'y')
    hold off
    % Get frame and save in video
    F = getframe;
    if i == 1
        xlen = size(F.cdata, 1);
        ylen = size(F.cdata, 2);
    else
        F.cdata = imresize(F.cdata, [xlen, ylen]);
    end
    writeVideo(writerObj, F);
end
close(writerObj)


