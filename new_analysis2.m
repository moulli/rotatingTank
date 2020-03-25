clear; close all; clc
addpath(genpath('/home/ljp/Science/balancoire_data/balancoire_matlab_files'))


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
numframes = 1000; %floor(vid.Duration * vid.FrameRate);
width = vid.Width;
height = vid.Height;


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
[x, y] = meshgrid(-(center(2)-1):(width-center(2)), -(center(1)-1):(height-center(1))); % meshgrid for mask
% ATTENTION: there is a confusion between center(1) and center(2) in the
% above equation. Lower, when defining another mask, the y-center was in
% the left part of meshgrid, and the x-center was in the right part of the
% meshgrid, so I switched it as well here to be sure...

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
fpts = cell(numframes, 1);
maskused = zeros(numframes, 1);
for frame = 1:numframes
    im = readFrame(vid);
    if frame == 1
        first_mask = 10; % arbitrary number for first mask
        maskused(frame) = first_mask;
        new_im = (background - im) .* cmasks(:, :, first_mask);
        [dpframex, dpframey] = find(new_im == max(new_im(:)));
        dpframex = dpframex(1);
        dpframey = dpframey(1);
        dp(frame, :) = [dpframex, dpframey];
    else
        % Chose appropriate mask
        nmask = size(cmasks, 3);
        dpmaskx = zeros(nmask, 1);
        dpmasky = zeros(nmask, 1);
        pixmask = zeros(nmask, 1);
        for mask = 1:nmask      
            % Get highest pixel
            new_im = (background - im) .* cmasks(:, :, mask);
            [dpframex, dpframey] = find(new_im == max(new_im(:)));
            dpmaskx(mask) = dpframex(1);
            dpmasky(mask) = dpframey(1);
            pixmask(mask) = max(new_im(:));
        end
        % Take mask with highest pixel value, and displacement not above threshold
        [~, pixind] = sort(pixmask, 'descend');
        for pixi = 1:length(pixind)
            pix = pixind(pixi);
            diff_disp = sqrt((dpmaskx(pix)-dp(frame-1, 1)).^2 + (dpmasky(pix)-dp(frame-1, 2)).^2);
            if abs(diff_disp) < max_diff_disp
                dpframex = dpmaskx(pix);
                dpframey = dpmasky(pix);
                dp(frame, :) = [dpframex, dpframey];
                ok_pix = min([pix+3, max(pixind)]);
                new_im = (background - im) .* cmasks(:, :, ok_pix);
                maskused(frame) = ok_pix;
                break
            end
            if pixi == length(pixind)
                pix = pixind(1);
                dpframex = dpmaskx(pix);
                dpframey = dpmasky(pix);
                dp(frame, :) = [dpframex, dpframey];
                new_im = (background - im) .* cmasks(:, :, pix);
                maskused(frame) = pix;
            end
        end
    end
    % Then make a mask for neighbourhood
    rnei = 45;
    [xnei, ynei] = meshgrid(-(dpframey-1):(width-dpframey), -(dpframex-1):(height-dpframex)); 
    masknei = uint8((xnei.^2 + ynei.^2) <= rnei^2);
    % Identifying important points
    new_im = masknei .* new_im;
    num_segments = 5;
    tail_length = 30;
    body_length = 0;
    inertia = 0.5;
    num_pix1 = 10;
    num_pix2 = 150;
    initial_box = 0.6;
    new_im = (255 - new_im);
    segment_pts = segmentTracking_WM(new_im, 'num_segments', num_segments, ...
                        'tail_length', tail_length, 'body_length', body_length, ...
                        'inertia', inertia, 'num_pix1', num_pix1, ...
                        'num_pix2', num_pix2, 'initial_box', initial_box);
    fpts{frame} = segment_pts;
end

% Plot fish journey
figure
hold on
image(im, 'CDataMapping', 'scaled')
plot(dp(:, 2), dp(:, 1))
axis equal


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
    % Treat image
    new_im = (background - im) .* cmasks(:, :, maskused(frame));
    [xnei, ynei] = meshgrid(-(dp(frame, 2)-1):(width-dp(frame, 2)), -(dp(frame, 1)-1):(height-dp(frame, 1))); 
    masknei = uint8((xnei.^2 + ynei.^2) <= rnei^2);
    new_im = 255 - masknei .* new_im;
    % Plot
    hold on
    image(new_im, 'CDataMapping', 'scaled')
    plot(ptss(2:2:end, frame), ptss(1:2:end, frame), 'r')
    axis equal
    pause(0.01)
    hold off
end

figure
subplot(2, 1, 1)
hold on
plot(ptss(1, :))
plot(ptss(end-1, :))
subplot(2, 1, 2)
hold on
plot(ptss(2, :))
plot(ptss(end, :))

figure
hold on
plot(diff(ptss(end-1, :)-ptss(1, :)))
plot(diff(ptss(end, :)-ptss(2, :)))

figure
ptss1 = sqrt(ptss(1, :).^2 + ptss(2, :).^2);
ptss2 = sqrt(ptss(end-1, :).^2 + ptss(end, :).^2);
plot(diff(ptss2 - ptss1))


%% How to better get fish image

% imbinarize
btemp = imbinarize(new_im);
figure
imshow(btemp)

% Fit exponential
temp = double(new_im);
temp(temp == 255) = [];
x = min(temp):254;
y = zeros(size(x));
for i = 1:length(x)
    y(i) = sum(temp == x(i));
end
efit = fit(x', y', 'exp1');
figure
plot(efit, x, y)
yg = efit.a * exp(efit.b .* x);
for i = length(yg):-1:1
    if yg(i) < 10
        thresh = i;
        break
    end
end
figure
new_im = double(new_im);
bin_im = (new_im .* (new_im <= x(thresh)));
bin_im(bin_im == 0) = 255;
image(bin_im, 'CDataMapping', 'scaled')
axis equal






