clear; close all; clc
addpath(genpath('/home/ljp/Science/balancoire_data/balancoire_matlab_files'))
addpath(genpath('/home/ljp/Science/Hippolyte/rotatingTank'))
folderpath = '/home/ljp/Science/rotatingTank_vidMatrices';


%% Load structure of interest

% Specify parameters
fish = 17;
light = 'On';
rot = 4;

% Load washingOut structure
path = fullfile(folderpath, strcat('washingOut_', num2str(fish), '_', light, '_', num2str(rot), '.mat'));
load(path, 'washingOut');


%% Plot fish position

% Animation
figure
inst = 1;
for frame = inst:washingOut.parameters.numframes
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
dt = 1 / washingOut.parameters.framerate;
time = (0:1:(washingOut.parameters.numframes-1)) * dt;
subplot(2, 1, 1)
hold on
for i = 1:washingOut.tracking.parameters.num_segments+1
    tseriex = reshape(washingOut.tracking.ptss(i, 1, :), [washingOut.parameters.numframes, 1]);
    plot(time, tseriex)
end
% plot(time, washingOut.tracking.modified*200)
title('x-axis displacement', 'Interpreter', 'latex')
xlabel('time [s]', 'Interpreter', 'latex')
ylabel('x-coordinate', 'Interpreter', 'latex')
subplot(2, 1, 2)
hold on
for i = 1:washingOut.tracking.parameters.num_segments+1
    tseriey = reshape(washingOut.tracking.ptss(i, 2, :), [washingOut.parameters.numframes, 1]);
    plot(time, tseriey)
end
% plot(time, washingOut.tracking.modified*200)
title('y-axis displacement', 'Interpreter', 'latex')
xlabel('time [s]', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')


%% Plot displacements differential

figure
dt = 1 / washingOut.parameters.framerate;
time = (0:1:(washingOut.parameters.numframes-1)) * dt;
subplot(2, 1, 1)
hold on
for i = 1:washingOut.tracking.parameters.num_segments+1
    tseriex = reshape(washingOut.tracking.ptss(i, 1, :), [washingOut.parameters.numframes, 1]);
    plot(time(1:end-1), diff(tseriex))
end
% plot(time, washingOut.tracking.modified*200)
title('x-axis displacement', 'Interpreter', 'latex')
xlabel('time [s]', 'Interpreter', 'latex')
ylabel('x-coordinate', 'Interpreter', 'latex')
subplot(2, 1, 2)
hold on
for i = 1:washingOut.tracking.parameters.num_segments+1
    tseriey = reshape(washingOut.tracking.ptss(i, 2, :), [washingOut.parameters.numframes, 1]);
    plot(time(1:end-1), diff(tseriey))
end
% plot(time, washingOut.tracking.modified*200)
title('y-axis displacement', 'Interpreter', 'latex')
xlabel('time [s]', 'Interpreter', 'latex')
ylabel('y-coordinate', 'Interpreter', 'latex')


%% Plot angle

figure
Angles = zeros(washingOut.parameters.numframes, washingOut.tracking.parameters.num_segments);
for frame = 1:washingOut.parameters.numframes
    [angles, angle0] = getTrackingAngles(washingOut.tracking.ptss(:, :, frame));
    Angles(frame, :) = [angle0, angles'];
end
total_angle = sum(Angles, 2);
plot(time, total_angle)
title('Sum of segments angles', 'Interpreter', 'latex')
xlabel('time [s]', 'Interpreter', 'latex')
ylabel('total angle', 'Interpreter', 'latex')


%% Plot bouts detector

figure
tseriex = reshape(washingOut.tracking.ptss(1, 1, :), [washingOut.parameters.numframes, 1]);
dtseriex = diff(tseriex);
tseriey = reshape(washingOut.tracking.ptss(1, 2, :), [washingOut.parameters.numframes, 1]);
dtseriey = diff(tseriey);
tserie = sqrt(dtseriex.^2 + dtseriey.^2);
hold on
plot(time(1:end-1), tserie)
plot(time(1:end-1), tserie >= 4)
title('Sum of segments angles', 'Interpreter', 'latex')
xlabel('time [s]', 'Interpreter', 'latex')
ylabel('total angle', 'Interpreter', 'latex')



