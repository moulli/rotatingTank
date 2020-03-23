clear; close all; clc


%% Get video

vidpath = '/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish1/Light On/Films/Fish On_1_12.avi';
vidpath = '/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish13/Light Off/Films/Fish Off_13_6.avi';
vidpath = '/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish18/Light On/Films/Fish On_18_0.avi';
vid = VideoReader(vidpath);
numframes = floor(vid.Duration * vid.FrameRate);
width = vid.Width;
height = vid.Height;


%% Get circle of tank

% Redefine video reader
vid = VideoReader(vidpath);

% Initialize background image
background = 0; 

% Add frame every dframe
dframe = 10;
for frame = 1:dframe:numframes
    im = read(vid, frame);
    background = background + double(im(:, :, 1));
end

% Compute final background
background = uint8(background./(numframes/dframe));

% Plot circles found
[centers, radii] = imfindcircles(background, [round((height)/2.3), round((height)/2)], 'ObjectPolarity', 'dark', 'Sensitivity', 0.99);
radius = max(radii);
center = centers(radius == radii, :);
% Arbitrary coefficient to get rid of tank borders
radius_coefficient = 0.95;
radius = radius_coefficient * radius;

% Plot results
figure
image(background, 'CDataMapping', 'scaled')
viscircles(center, radius, 'EdgeColor', 'b');
axis equal

% Create mask
[x, y] = meshgrid(-(center(1)-1):(width-center(1)), -(center(2)-1):(height-center(2)));
cmask = uint8((x.^2 + y.^2) <= radius^2);
figure
image(im .* cmask, 'CDataMapping', 'scaled')
axis equal


%% Get lighter point at each frame using background and mask

% Redefine video reader
vid = VideoReader(vidpath);

% Loop over all frames
dp = zeros(numframes, 2);
for frame = 1:numframes
    im = readFrame(vid);
    new_im = (background - im) .* cmask;
    [dpframex, dpframey] = find(new_im == max(new_im(:)));
    dp(frame, :) = [dpframex(1), dpframey(1)];
end

% Plot fish journey
figure
hold on
image(im, 'CDataMapping', 'scaled')
plot(dp(:, 2), dp(:, 1))
axis equal
% scatter(dpind(:, 1), dpind(:, 2))


%% There is a problem with the contour
% Either it is too large, and the algorithm sometimes locates the fish in
% this area, where it should not be, resulting in position jumps, either it
% is too short, and the algorithm fails to locate the fish when it it
% swimming in this area. Let's work on an adaptive contour radius to solve
% this problem.

% Getting largest radius and associated center
radius = max(radii);
center = centers(radius == radii, :);

% Getting averaged pixel value from center of video
[x, y] = meshgrid(-(center(1)-1):(width-center(1)), -(center(2)-1):(height-center(2)));
cmaskmiddle = uint8((x.^2 + y.^2) <= (height/5)^2);
pixvalue = sum(sum(background .* cmaskmiddle)) / sum(sum(cmaskmiddle));

% Computing averaged pixel value for different radii
mult_coeffs = 0.80:0.001:1.2;
pixradius = zeros(size(mult_coeffs));
for cmul = 1:length(mult_coeffs)
    cmasktemps = uint8((x.^2 + y.^2) <= (mult_coeffs(cmul)*radius)^2);
    pixradius(cmul) = sum(sum(background .* cmasktemps)) / sum(sum(cmasktemps));
end
% Plot pixradius and difference of pixradius
figure
subplot(2, 1, 1)
plot(mult_coeffs, pixradius)
title('Evolution of averaged pixel value depending on mask radius', 'Interpreter', 'latex')
subplot(2, 1, 2)
plot(mult_coeffs(1:end-1), diff(pixradius), 'r')
title('Differential of averaged pixel value depending on mask radius', 'Interpreter', 'latex')

% Get radius before minimum differential of pixradius
[~, mindiff] = min(diff(pixradius));
radius_coefficient = mult_coeffs(mindiff) - 0.03;
radius = radius_coefficient * radius;

% Plot results
figure
image(background, 'CDataMapping', 'scaled')
viscircles(center, radius, 'EdgeColor', 'b');
axis equal

% Create mask
[x, y] = meshgrid(-(center(1)-1):(width-center(1)), -(center(2)-1):(height-center(2)));
cmask = uint8((x.^2 + y.^2) <= radius^2);


%% Solution: store different masks, and update it if weird fish movement

% Redefine video reader
vid = VideoReader(vidpath);

% Initialize background image
background = 0; 

% Add frame every dframe
dframe = 10;
for frame = 1:dframe:numframes
    im = read(vid, frame);
    background = background + double(im(:, :, 1));
end

% Compute final background
background = uint8(background./(numframes/dframe));

% Get circle with highest radius
[centers, radii] = imfindcircles(background, [round((height)/2.3), round((height)/2)], 'ObjectPolarity', 'dark', 'Sensitivity', 0.99);
radius = max(radii);
center = centers(radius == radii, :);

% Computing averaged pixel value for different radii
mult_coeffs = 0.80:0.001:1.2;
pixradius = zeros(size(mult_coeffs));
for cmul = 1:length(mult_coeffs)
    cmasktemps = uint8((x.^2 + y.^2) <= (mult_coeffs(cmul)*radius)^2);
    pixradius(cmul) = sum(sum(background .* cmasktemps)) / sum(sum(cmasktemps));
end

% Get radius before minimum differential of pixradius
[~, mindiff] = min(diff(pixradius));
radius_coeff = mult_coeffs(mindiff);

% Create masks
radius_coeff_inc = 0.01;
radii = ((radius_coeff-10*radius_coeff_inc):radius_coeff_inc:(radius_coeff+5*radius_coeff_inc)) .* radius;
[x, y] = meshgrid(-(center(1)-1):(width-center(1)), -(center(2)-1):(height-center(2)));
cmasks = zeros(height, width, length(radii));
for r = 1:length(radii)
    cmaskr = ((x.^2 + y.^2) <= radii(r)^2);
    cmasks(:, :, r) = cmaskr;
end
cmasks = uint8(cmasks);


% Redefine video reader
vid = VideoReader(vidpath);

% Define maximum displacement differential to switch mask in needed
max_diff_disp = 20;

% Loop over all frames
dp = zeros(numframes, 2);
for frame = 1:numframes
    im = readFrame(vid);
    if frame == 1
        first_mask = 8; % arbitrary number for first mask
        new_im = (background - im) .* cmasks(:, :, first_mask);
        [dpframex, dpframey] = find(new_im == max(new_im(:)));
        dp(frame, :) = [dpframex(1), dpframey(1)];
    else
        for mask = 1:size(cmasks, 3)       
            % Get highest pixel
            new_im = (background - im) .* cmasks(:, :, mask);
            [dpframex, dpframey] = find(new_im == max(new_im(:)));
            % Check displacement is not above threshold
            dpframex = dpframex(1);
            dpframey = dpframey(1);
            diff_disp = sqrt((dpframex-dp(frame-1, 1)).^2 + (dpframey-dp(frame-1, 2)).^2);
            if abs(diff_disp) < max_diff_disp || mask == size(cmasks, 3)
                dp(frame, :) = [dpframex, dpframey];
                break
            end
        end
    end
end

% Plot fish journey
figure
hold on
image(im, 'CDataMapping', 'scaled')
plot(dp(:, 2), dp(:, 1))
axis equal






