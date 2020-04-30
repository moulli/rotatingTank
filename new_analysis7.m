clear; close all; clc
addpath(genpath('/home/ljp/Science/balancoire_data/balancoire_matlab_files'))
addpath(genpath('/home/ljp/Science/Hippolyte/rotatingTank'))
folderpath = '/home/ljp/Science/Hippolyte/rotatingTank_vidMatrices';


%% Loop over all files for analysis

fish = 1:18;
light = {'Off', 'On'};
rot = 0:2:12;

savecell = cell(length(fish), length(light), length(rot));

for f = fish
    for l = light
        for r = rot
            % Load washingOut
            matname = strcat('washingOut_', num2str(f), '_', l{1}, '_', num2str(r), '.mat');
            path = fullfile(folderpath, matname);
            load(path, 'washingOut');
            % Define statistics structure
            statistics = struct;
            % First distance made at each frame
            statistics.distance = sqrt(sum((washingOut.images.dp(2:end, :) - washingOut.images.dp(1:end-1, :)).^2, 2));
            % Then distance from center of tank and angle in tank
            dp_centered = washingOut.images.dp - washingOut.preprocessing.center;
            [statistics.tank_angle, statistics.radius] = cart2pol(dp_centered(:, 1), dp_centered(:, 2)); 
            statistics.tank_angle = 180/pi * statistics.tank_angle;
            % Finally distance without rotation
            rad_per_frame = pi/180 * r / washingOut.parameters.framerate;
            rad_per_frame = rad_per_frame * (0:(washingOut.parameters.numframes-1))';
            [dp_rot_adjust(:, 1), dp_rot_adjust(:, 2)] = pol2cart(statistics.tank_angle*pi/180 - rad_per_frame, statistics.radius);
            statistics.distance_notrot = sqrt(sum((dp_rot_adjust(2:end, :) - dp_rot_adjust(1:end-1, :)).^2, 2));
            
            if isequal(l{1}, 'Off')
                ln = 1;
            else
                ln = 2;
            end
            savecell{f, ln, find(rot == r)} = statistics;
            disp([f, ln, r]);
        end
    end
end


%% Plotting boxplots

dist = zeros(length(fish), length(rot), length(light));
tang = zeros(length(fish), length(rot), length(light));
radi = zeros(length(fish), length(rot), length(light));

for f = 1:length(fish)
    for l = 1:length(light)
        for r = 1:length(rot)
            dist(f, r, l) = mean(savecell{f, l, r}.distance);
            tang(f, r, l) = mean(savecell{f, l, r}.tank_angle);
            radi(f, r, l) = mean(savecell{f, l, r}.radius);
        end
    end
end

indexes = ones(length(fish), 1) .* (1:length(rot));
indexes = indexes(:);

figure
hold on
dist_off = dist(:, :, 1);
dist_on = dist(:, :, 2);
boxplot(dist_off(:), indexes, 'Colors', 'k')
boxplot(dist_on(:), indexes, 'Colors', 'r')
axis([0.5, 7.5, 0, 17])

figure
hold on
tang_off = tang(:, :, 1);
tang_on = tang(:, :, 2);
boxplot(tang_off(:), indexes, 'Colors', 'k')
boxplot(tang_on(:), indexes, 'Colors', 'r')

figure
hold on
radi_off = radi(:, :, 1);
radi_on = radi(:, :, 2);
boxplot(radi_off(:), indexes, 'Colors', 'k')
boxplot(radi_on(:), indexes, 'Colors', 'r')





            