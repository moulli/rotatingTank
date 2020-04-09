clear; close all; clc
addpath(genpath('/home/ljp/Science/balancoire_data/balancoire_matlab_files'))
addpath(genpath('/home/ljp/Science/Hippolyte/rotatingTank'))
folderpath = '/home/ljp/Science/rotatingTank_vidMatrices';

foldertemp = '/home/ljp/Science/rotatingTank_vidMatrices/matrices_old';
temp = dir(foldertemp);
for file = 3:length(temp)
    % load
    pathload = fullfile(foldertemp, temp(file).name);
    load(pathload, 'washingOut')
    % modify structure
    washingOut2 = struct;
    washingOut2.parameters = washingOut.parameters;
    washingOut2.parameters.duration = washingOut.duration;
    washingOut2.parameters.framerate = washingOut.framerate;
    washingOut2.preprocessing = washingOut.preprocessing;
    washingOut2.images = washingOut.images;
    washingOut2.tracking = washingOut.tracking;
    washingOut = washingOut2;
    % add center for preprocessing
    m = washingOut.preprocessing.masks(:, :, 1);
    i = find(m == 1);
    [x, y] = ind2sub(size(m), i);
    xc = (min(x)+max(x))/2;
    yc = (min(y)+max(y))/2;
    washingOut.preprocessing.center = [xc, yc];
    % save new structure
    pathsave = fullfile(folderpath, temp(file).name);
    save(pathsave, 'washingOut');
    % print
    disp(file)
end