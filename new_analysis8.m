clear; close all; clc
addpath('/home/ljp/Science/Loann/Loann Codes and Content/Scripts Matlab')



%% Save all required information from Loann matrices

% Define parameters

fish = 3:18;
light = {' Off', ' On'};
rot = 0:2:12;


% Define output structure

MegaFish = cell(length(fish), length(light), length(rot));


% Loop over all fish, light and rotation

megapath = '/home/ljp/Science/Hippolyte/MegaFish.mat';
if exist(megapath, 'file')
    load(megapath);
else

    warning('off')
    for floop = fish
        for lloop = light
            for rloop = rot

                %% Get path to matrix

                matpath = strcat('/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish', ...
                                  num2str(floop), '/Light', lloop{1}, '/Matlab Data/Fish FILTERED', lloop{1}, '_', ...
                                  num2str(floop), '_', num2str(rloop), '.mat');


                %% Load data

                load(matpath);


                %% Create video reader

                vidpath = strcat('/home/ljp/Science/Loann/Results Loann/Fish Data light on and off/Fish', ...
                                  num2str(floop), '/Light', lloop{1}, '/Films/Fish' , lloop{1}, '_', ...
                                  num2str(floop), '_', num2str(rloop), '.avi');
                v = VideoReader(vidpath);


                %% Get info in structure

                Fish = Fish_Computation(Angle, speed, angle_r, angle_l, t, v, posx, posy);


                %% Put structure in mega structure

                if isequal(lloop{1}, ' Off')
                    ln = 1;
                else
                    ln = 2;
                end
                MegaFish{(fish==floop), ln, (rot==rloop)} = Fish;


            end
        end
    end
    warning('on')
    
end



%%

floop = randperm(length(fish), 1);
lloop = randperm(length(light), 1);
rloop = randperm(length(rot), 1);
Fish = MegaFish{floop, lloop, rloop};
figure
hold on
plot(Fish.angle)
scatter(Fish.Bouts_Starts, Fish.angle(Fish.Bouts_Starts))
scatter(Fish.Bouts_Ends, Fish.angle(Fish.Bouts_Ends))
axis([1, length(Fish.angle), -180, 180])
title('angle of the fish with bouts starts (orange) and bouts ends (yellow)', 'Interpreter', 'latex')
xlabel('frame', 'Interpreter', 'latex')
ylabel('angle [degree]', 'Interpreter', 'latex')



%%

info_time = zeros(0, 7);
info_bout = cell(0, 2);
info_video = cell(0, 3);
for floop = 1:length(fish)
    for lloop = 1:length(light)
        for rloop = 1:length(rot)
            angle = MegaFish{floop, lloop, rloop}.angle;
            sbouts = MegaFish{floop, lloop, rloop}.Bouts_Starts;
            ebouts = MegaFish{floop, lloop, rloop}.Bouts_Ends;
            ind = 1;
            while ind <= length(angle)
                % if nan get to next index
                if isnan(angle(ind))
                    ind = ind + 1;
                    continue
                end
                % else count to next nan
                chaine = 1;
                keepgoing = true;
                while keepgoing
                    % stop if end of angle vector or next value is nan
                    if ind + chaine > length(angle) || isnan(angle(ind+chaine))
                        keepgoing = false;
                    % else increment chaine
                    else
                        chaine = chaine + 1;
                    end
                end
                num_bouts = sum(ismember(ind:(ind+chaine-1), sbouts));
                addinfo = [floop, lloop, rloop, chaine, num_bouts, ind, ind+chaine-1];
                info_time = cat(1, info_time, addinfo);
                keepbout = ismember(sbouts, ind:(ind+chaine-1));
                info_bout = [info_bout; [{sbouts(keepbout)}, {ebouts(keepbout)}]];
                info_video = [info_video; [{fish(floop)}, {light{lloop}}, {rot(rloop)}]];
                ind = ind+chaine;
            end
        end
    end
end

% Define structure
info_fish = struct;
info_fish.info_time = info_time;
info_fish.info_bout = info_bout;
info_fish.info_video = info_video;

% Plot lengths of sequences and number of bouts by quantile
qs = 0:0.01:1;  
numq = zeros(length(qs), 1);
numb = zeros(length(qs), 1);
for q = 1:length(qs)
    numq(q) = quantile(info_fish.info_time(:, 4), qs(q));
    numb(q) = quantile(info_fish.info_time(:, 5), qs(q));
end  
figure
subplot(2, 1, 1)
plot(qs, numq)
axis([0, 1, 0, 100])
title('lenght of sequence without nan value, depending on quantile for number of sequences', 'Interpreter', 'latex')
xlabel('quantile', 'Interpreter', 'latex')
ylabel('length of sequence', 'Interpreter', 'latex')
subplot(2, 1, 2)
plot(qs, numb)
axis([0, 1, 0, 10])
title('number of bouts in sequence without nan value, depending on quantile for number of sequences', 'Interpreter', 'latex')
xlabel('quantile', 'Interpreter', 'latex')
ylabel('number of bouts in sequence', 'Interpreter', 'latex')


%% Multilinear regression using Schopik's paper

% First isolate time series with at least one bout
oneb = (info_fish.info_time(:, 5) > 0);
onetime = info_fish.info_time(oneb, :);
onebout = info_fish.info_bout(oneb, :);
onelen = sum(oneb);

% - angle after bout:
outangle = zeros(0, 1);
% - angle correction:
corangle = zeros(0, 1);
% - absolute angle before bout:
absangle = zeros(0, 1);
% - absolute angle differential before bout:
absdiffang = zeros(0, 1);
% - angle differential before bout:
diffang = zeros(0, 1);

% Then get for each bout all the important information
for o = 1:onelen
    time = onetime(o, :);
    angle = MegaFish{time(1), time(2), time(3)}.angle;
    startpts = onebout{o, 1};
    endpts = onebout{o, 2};
    for a = 1:length(endpts)
        outangle = cat(1, outangle, angle(endpts(a)));
        corangle = cat(1, corangle, angle(endpts(a))-angle(startpts(a)));
        absangle = cat(1, absangle, abs(angle(startpts(a))));
        absdiffang = cat(1, absdiffang, abs(angle(startpts(a))-angle(startpts(a)-1)));
        diffang = cat(1, diffang, angle(startpts(a))-angle(startpts(a)-1));
    end
end
boutlen = length(outangle);

% build regressors and outputs
y1 = (outangle - mean(outangle)) ./ std(outangle);
y2 = (corangle - mean(corangle)) ./ std(corangle);
X = zeros(boutlen, 4);
X(:, 1) = 1;
X(:, 2) = (absangle - mean(absangle)) ./ std(absangle);
X(:, 3) = (absdiffang - mean(absdiffang)) ./ std(absdiffang);
X(:, 4) = (diffang - mean(diffang)) ./ std(diffang);

% compute regressions
coeffs1 = (X' * X) \ X' * y1;
resid1 = y1 - X * coeffs1;
coeffs2 = (X' * X) \ X' * y2;
resid2 = y2 - X * coeffs2;


%% Multilinear regression using former bout information

% First isolate time series with at least one bout
oneb = (info_fish.info_time(:, 5) > 1);
onetime = info_fish.info_time(oneb, :);
onebout = info_fish.info_bout(oneb, :);
onelen = sum(oneb);

% - angle after bout:
outangle_l = zeros(0, 1);
outangle_r = zeros(0, 1);
% - angle correction:
corangle_l = zeros(0, 1);
corangle_r = zeros(0, 1);
% - absolute angle before bout:
absangle_l = zeros(0, 1);
absangle_r = zeros(0, 1);
% - absolute angle differential before bout:
absdiffang_l = zeros(0, 1);
absdiffang_r = zeros(0, 1);
% - angle differential before bout:
diffang_l = zeros(0, 1);
diffang_r = zeros(0, 1);

% Then get for each bout all the important information
for o = 1:onelen
    time = onetime(o, :);
    angle = MegaFish{time(1), time(2), time(3)}.angle;
    startpts = onebout{o, 1};
    endpts = onebout{o, 2};
    for a = 2:length(endpts)
        if angle(startpts(a)) < 0
            outangle_l = cat(1, outangle_l, angle(endpts(a)));
            corangle_l = cat(1, corangle_l, angle(endpts(a))-angle(startpts(a)));
            absangle_l = cat(1, absangle_l, abs(angle(startpts(a))));
            absdiffang_l = cat(1, absdiffang_l, abs(angle(startpts(a))-angle(startpts(a)-1)));
            diffang_l = cat(1, diffang_l, angle(startpts(a))-angle(startpts(a)-1));
        else
            outangle_r = cat(1, outangle_r, angle(endpts(a)));
            corangle_r = cat(1, corangle_r, angle(endpts(a))-angle(startpts(a)));
            absangle_r = cat(1, absangle_r, abs(angle(startpts(a))));
            absdiffang_r = cat(1, absdiffang_r, abs(angle(startpts(a))-angle(startpts(a)-1)));
            diffang_r = cat(1, diffang_r, angle(startpts(a))-angle(startpts(a)-1));
        end
    end
end
boutlen_l = length(outangle_l);
boutlen_r = length(outangle_r);

% build regressors and outputs for left
y1_l = (outangle_l - mean(outangle_l)) ./ std(outangle_l);
y2_l = (corangle_l - mean(corangle_l)) ./ std(corangle_l);
X_l = zeros(boutlen_l, 4);
X_l(:, 1) = 1;
X_l(:, 2) = (absangle_l - mean(absangle_l)) ./ std(absangle_l);
X_l(:, 3) = (absdiffang_l - mean(absdiffang_l)) ./ std(absdiffang_l);
X_l(:, 4) = (diffang_l - mean(diffang_l)) ./ std(diffang_l);
% build regressors and outputs for right
y1_r = (outangle_r - mean(outangle_r)) ./ std(outangle_r);
y2_r = (corangle_r - mean(corangle_r)) ./ std(corangle_r);
X_r = zeros(boutlen_r, 4);
X_r(:, 1) = 1;
X_r(:, 2) = (absangle_r - mean(absangle_r)) ./ std(absangle_r);
X_r(:, 3) = (absdiffang_r - mean(absdiffang_r)) ./ std(absdiffang_r);
X_r(:, 4) = (diffang_r - mean(diffang_r)) ./ std(diffang_r);

% compute regressions for left
coeffs1_l = (X_l' * X_l) \ X_l' * y1_l;
resid1_l = y1_l - X_l * coeffs1_l;
coeffs2_l = (X_l' * X_l) \ X_l' * y2_l;
resid2_l = y2_l - X_l * coeffs2_l;
% compute regressions for right
coeffs1_r = (X_r' * X_r) \ X_r' * y1_r;
resid1_r = y1_r - X_r * coeffs1_r;
coeffs2_r = (X_r' * X_r) \ X_r' * y2_r;
resid2_r = y2_r - X_r * coeffs2_r;

% NOTA BENE:
% instead of splitting dataset between left and right, we could have simply
% added an orientation parameter, equal to -1 if fish is facing left before
% bout, and 1 if fish is facing right. The main problem of our former
% algorithm (which used absolute angle as in Schoppik paper), was that no
% linear regression could accurately provide orientation post bout due to
% lack of information on orientation pre bout.

figure
subplot(2, 1, 1)
hold on
plot(y1_l(1:250))
plot(X_l(1:250, :) * coeffs1_l)
title('Left facing fish', 'Interpreter', 'latex')
xlabel('example', 'Interpreter', 'latex')
ylabel('normalized post-bout angle', 'Interpreter', 'latex')
legend('actual data', 'regression data')
subplot(2, 1, 2)
hold on
plot(y1_r(1:250))
plot(X_r(1:250, :) * coeffs1_r)
title('Right facing fish', 'Interpreter', 'latex')
legend('actual data', 'regression data')
xlabel('example', 'Interpreter', 'latex')
ylabel('normalized post-bout angle', 'Interpreter', 'latex')


%% Comparing different models

% First isolate time series with at least one bout
oneb = (info_fish.info_time(:, 5) > 1);
onetime = info_fish.info_time(oneb, :);
onebout = info_fish.info_bout(oneb, :);
onelen = sum(oneb);

% - angle after bout:
outangle_l = zeros(0, 1);
outangle_r = zeros(0, 1);
% - absolute angle before bout:
absangle_l = zeros(0, 1);
absangle_r = zeros(0, 1);
% - absolute angle differential before bout:
absdiffang_l = zeros(0, 1);
absdiffang_r = zeros(0, 1);
% - angle differential before bout:
diffang_l = zeros(0, 1);
diffang_r = zeros(0, 1);
% - angle after bout on former bout
formout_l = zeros(0, 1);
formout_r = zeros(0, 1);

% Then get for each bout all the important information
for o = 1:onelen
    time = onetime(o, :);
    angle = MegaFish{time(1), time(2), time(3)}.angle;
    startpts = onebout{o, 1};
    endpts = onebout{o, 2};
    for a = 2:length(endpts)
        if angle(startpts(a)) < 0
            outangle_l = cat(1, outangle_l, angle(endpts(a)));
            absangle_l = cat(1, absangle_l, abs(angle(startpts(a))));
            absdiffang_l = cat(1, absdiffang_l, abs(angle(startpts(a))-angle(startpts(a)-1)));
            diffang_l = cat(1, diffang_l, angle(startpts(a))-angle(startpts(a)-1));
            formout_l = cat(1, formout_l, angle(endpts(a-1)));
        else
            outangle_r = cat(1, outangle_r, angle(endpts(a)));
            absangle_r = cat(1, absangle_r, abs(angle(startpts(a))));
            absdiffang_r = cat(1, absdiffang_r, abs(angle(startpts(a))-angle(startpts(a)-1)));
            diffang_r = cat(1, diffang_r, angle(startpts(a))-angle(startpts(a)-1));
            formout_r = cat(1, formout_r, angle(endpts(a-1)));
        end
    end
end
boutlen_l = length(outangle_l);
boutlen_r = length(outangle_r);

test_num = 500;
% create training and test sets for left
mix_l = randperm(boutlen_l);
train_l = mix_l(1:boutlen_l-test_num);
test_l = mix_l(boutlen_l-test_num+1:end);
% create training and test sets for right
mix_r = randperm(boutlen_r);
train_r = mix_r(1:boutlen_r-test_num);
test_r = mix_r(boutlen_r-test_num+1:end);

% build regressors and outputs for left
y_l = outangle_l(train_l);
Xall_l = zeros(boutlen_l, 4);
Xall_l(:, 1) = 1;
Xall_l(:, 2) = (absangle_l - mean(absangle_l)) ./ std(absangle_l);
Xall_l(:, 3) = (absdiffang_l - mean(absdiffang_l)) ./ std(absdiffang_l);
Xall_l(:, 4) = (diffang_l - mean(diffang_l)) ./ std(diffang_l);
X_l = Xall_l(train_l, :);
% build regressors and outputs for right
y_r = outangle_r(train_r);
Xall_r = zeros(boutlen_r, 4);
Xall_r(:, 1) = 1;
Xall_r(:, 2) = (absangle_r - mean(absangle_r)) ./ std(absangle_r);
Xall_r(:, 3) = (absdiffang_r - mean(absdiffang_r)) ./ std(absdiffang_r);
Xall_r(:, 4) = (diffang_r - mean(diffang_r)) ./ std(diffang_r);
X_r = Xall_r(train_r, :);

% compute regressions for left
coeffs_l = (X_l' * X_l) \ X_l' * y_l;
resid_l = outangle_l(test_l) - Xall_l(test_l, :) * coeffs_l;
% compute regressions for right
coeffs_r = (X_r' * X_r) \ X_r' * y_r;
resid_r = outangle_r(test_r) - Xall_r(test_r, :) * coeffs_r;

figure
subplot(2, 1, 1)
hold on
plot(outangle_l(test_l))
plot(Xall_l(test_l, :) * coeffs_l)
title('Left facing fish regression on test set', 'Interpreter', 'latex')
xlabel('example', 'Interpreter', 'latex')
ylabel('post-bout angle', 'Interpreter', 'latex')
legend('actual data', 'regression data')
subplot(2, 1, 2)
hold on
plot(outangle_r(test_r))
plot(Xall_r(test_r, :) * coeffs_r)
title('Right facing fish regression on test set', 'Interpreter', 'latex')
legend('actual data', 'regression data')
xlabel('example', 'Interpreter', 'latex')
ylabel('post-bout angle', 'Interpreter', 'latex')


%% Simplifying using function

% First isolate time series with at least one bout
oneb = (info_fish.info_time(:, 5) > 1);
onetime = info_fish.info_time(oneb, :);
onebout = info_fish.info_bout(oneb, :);
onelen = sum(oneb);

% - angle after bout:
outangle_l = zeros(0, 1);
outangle_r = zeros(0, 1);
% - absolute angle before bout:
absangle_l = zeros(0, 1);
absangle_r = zeros(0, 1);
% - absolute angle differential before bout:
absdiffang_l = zeros(0, 1);
absdiffang_r = zeros(0, 1);
% - angle differential before bout:
diffang_l = zeros(0, 1);
diffang_r = zeros(0, 1);
% - angle after bout on former bout
formout_l = zeros(0, 1);
formout_r = zeros(0, 1);

% Then get for each bout all the important information
for o = 1:onelen
    time = onetime(o, :);
    angle = MegaFish{time(1), time(2), time(3)}.angle;
    startpts = onebout{o, 1};
    endpts = onebout{o, 2};
    for a = 2:length(endpts)
        if angle(startpts(a)) < 0 && -180 <= angle(startpts(a)) % -180 is because there is a weird point at -500 degrees
            outangle_l = cat(1, outangle_l, angle(endpts(a)));
            absangle_l = cat(1, absangle_l, abs(angle(startpts(a))));
            absdiffang_l = cat(1, absdiffang_l, abs(angle(startpts(a))-angle(startpts(a)-1)));
            diffang_l = cat(1, diffang_l, angle(startpts(a))-angle(startpts(a)-1));
            formout_l = cat(1, formout_l, angle(endpts(a-1)));
        elseif 0 <= angle(startpts(a)) && angle(startpts(a)) <= 180
            outangle_r = cat(1, outangle_r, angle(endpts(a)));
            absangle_r = cat(1, absangle_r, abs(angle(startpts(a))));
            absdiffang_r = cat(1, absdiffang_r, abs(angle(startpts(a))-angle(startpts(a)-1)));
            diffang_r = cat(1, diffang_r, angle(startpts(a))-angle(startpts(a)-1));
            formout_r = cat(1, formout_r, angle(endpts(a-1)));
        end
    end
end

test_prop = 0.1;
% regression for left using Schoppik's model
[train_coeffs_l, test_residuals_l, ~, ~, X_test_l, y_test_l, mix_l] = ...
            train_test_regress(test_prop, outangle_l, absangle_l, absdiffang_l, diffang_l);
[train_coeffs_r, test_residuals_r, ~, ~, X_test_r, y_test_r, mix_r] = ...
            train_test_regress(test_prop, outangle_r, absangle_r, absdiffang_r, diffang_r);
    
% adding former bouts end angles
[train_coeffs_l_f, test_residuals_l_f, ~, ~, X_test_l_f, y_test_l_f] = ...
            train_test_regress(test_prop, outangle_l, absangle_l, absdiffang_l, diffang_l, formout_l, mix_l);
[train_coeffs_r_f, test_residuals_r_f, ~, ~, X_test_r_f, y_test_r_f] = ...
            train_test_regress(test_prop, outangle_r, absangle_r, absdiffang_r, diffang_r, formout_r, mix_r);

% plotting        
figure
subplot(2, 1, 1)
hold on
plot(y_test_l)
plot(X_test_l * train_coeffs_l)
plot(y_test_l_f) % should be equal to y_test_l
plot(X_test_l_f * train_coeffs_l_f)
title('Left facing fish regression on test set', 'Interpreter', 'latex')
xlabel('example', 'Interpreter', 'latex')
ylabel('post-bout angle', 'Interpreter', 'latex')
legend('actual data', 'regression data', 'actual data', 'regression data with former angle')
subplot(2, 1, 2)
hold on
plot(y_test_r)
plot(X_test_r * train_coeffs_r)
plot(y_test_r_f) % should be equal to y_test_r
plot(X_test_r_f * train_coeffs_r_f)
title('Right facing fish regression on test set', 'Interpreter', 'latex')
legend('actual data', 'regression data', 'actual data', 'regression data with former angle')
xlabel('example', 'Interpreter', 'latex')
ylabel('post-bout angle', 'Interpreter', 'latex')


%% Computing F-tests for different regression configurations

test_prop = 0.1;

% define number of regressors and results cells
params = [1, 2, 2, 4, 5];
resids_l = cell(length(params), 1);
coeffs_l = cell(length(params), 1);
resids_r = cell(length(params), 1);
coeffs_r = cell(length(params), 1);
xtest_l = cell(length(params), 1);
xtest_r = cell(length(params), 1);

% fill these cells
[coeffs_l{1}, resids_l{1}, ~, ~, xtest_l{1}, yout_l, mix_l] = train_test_regress(test_prop, outangle_l);
[coeffs_r{1}, resids_r{1}, ~, ~, xtest_r{1}, yout_r, mix_r] = train_test_regress(test_prop, outangle_r);
[coeffs_l{2}, resids_l{2}, ~, ~, xtest_l{2}] = train_test_regress(test_prop, outangle_l, absangle_l, mix_l);
[coeffs_r{2}, resids_r{2}, ~, ~, xtest_r{2}] = train_test_regress(test_prop, outangle_r, absangle_r, mix_r);
[coeffs_l{3}, resids_l{3}, ~, ~, xtest_l{3}] = train_test_regress(test_prop, outangle_l, formout_l, mix_l);
[coeffs_r{3}, resids_r{3}, ~, ~, xtest_r{3}] = train_test_regress(test_prop, outangle_r, formout_r, mix_r);
[coeffs_l{4}, resids_l{4}, ~, ~, xtest_l{4}] = train_test_regress(test_prop, outangle_l, absangle_l, absdiffang_l, diffang_l, mix_l);
[coeffs_r{4}, resids_r{4}, ~, ~, xtest_r{4}] = train_test_regress(test_prop, outangle_r, absangle_r, absdiffang_r, diffang_r, mix_r);
[coeffs_l{5}, resids_l{5}, ~, ~, xtest_l{5}] = train_test_regress(test_prop, outangle_l, absangle_l, absdiffang_l, diffang_l, formout_l, mix_l);
[coeffs_r{5}, resids_r{5}, ~, ~, xtest_r{5}] = train_test_regress(test_prop, outangle_r, absangle_r, absdiffang_r, diffang_r, formout_r, mix_r);

% compute matrix of F-tests
ftests_l = zeros(length(params), length(params));
ftests_r = zeros(length(params), length(params));
for i = 1:length(params)
    for j = 1:length(params)
        if j > i
            fstat_l = ((sum(resids_l{i}.^2)-sum(resids_l{j}.^2))/(j-i))/ (sum(resids_l{j}.^2)/(length(yout_l)-j));
            ftests_l(i, j) = fpdf(fstat_l, j-i, length(yout_l)-j);
            fstat_r = ((sum(resids_r{i}.^2)-sum(resids_r{j}.^2))/(j-i))/ (sum(resids_r{j}.^2)/(length(yout_r)-j));
            ftests_r(i, j) = fpdf(fstat_r, j-i, length(yout_r)-j);
        end
    end
end
figure
subplot(1, 2, 1)
image(ftests_l, 'CDataMapping', 'scaled')
colorbar
axis equal
title('Hypothesis testing for regressor addition', 'Interpreter', 'latex')
subplot(1, 2, 2)
image(ftests_r, 'CDataMapping', 'scaled')
colorbar
axis equal
title('Hypothesis testing for regressor addition', 'Interpreter', 'latex')

% plot results
numpts = 100;
figure
subplot(2, 1, 1)
hold on
ytemp = yout_l(1:numpts);
[~, sort_ind] = sort(ytemp);
plot(yout_l(sort_ind), ':')
for i = 1:length(params)
    plot(xtest_l{i}(sort_ind, :) * coeffs_l{i})
end
title('Left facing fish regression on test set', 'Interpreter', 'latex')
legend('actual data', 'no parameters', 'with absolute angle', 'with former angle', 'schoppik model', 'schoppik model and former angle')
xlabel('example', 'Interpreter', 'latex')
ylabel('post-bout angle', 'Interpreter', 'latex')
subplot(2, 1, 2)
hold on
ytemp = yout_r(1:numpts);
[~, sort_ind] = sort(ytemp);
plot(yout_r(sort_ind), ':')
for i = 1:length(params)
    plot(xtest_r{i}(sort_ind, :) * coeffs_r{i})
end
title('Right facing fish regression on test set', 'Interpreter', 'latex')
legend('actual data', 'no parameters', 'with absolute angle', 'with former angle', 'schoppik model', 'schoppik model and former angle')
xlabel('example', 'Interpreter', 'latex')
ylabel('post-bout angle', 'Interpreter', 'latex')


