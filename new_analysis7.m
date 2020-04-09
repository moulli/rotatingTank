clear; close all; clc
addpath(genpath('/home/ljp/Science/balancoire_data/balancoire_matlab_files'))
addpath(genpath('/home/ljp/Science/Hippolyte/rotatingTank'))
folderpath = '/home/ljp/Science/rotatingTank_vidMatrices';


%% Loop over all files for analysis

fish = 1:18;
light = {'Off', 'On'};
rot = 0:2:12;

for f = fish
    for l = light
        for r = rot
            % Load washingOut
            matname = strcat('washingOut_', num2str(f), '_', l{1}, '_', num2str(r), '.mat');
            path = fullfile(folderpath, matname);
            load(path, 'washingOut');
            % Define statistics structure
            statistics = struct;
            % First distance made at each frame, with and without rotation, and distance to center of tank
            statistics.distance = sqrt(sum((washingOut.images.dp(2:end, :) - washingOut.images.dp(1:end-1, :)).^2, 2));
            dp_centered = washingOut.images.dp - washingOut.preprocessing.center;
            statistics.radius = sqrt(sum((dp_centered).^2, 2));
            statistics.tank_angle = atan(dp_centered(:, 2) ./ dp_centered(:, 1)) * 180 / pi;
        end
    end
end