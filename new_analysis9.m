clear; close all; clc
addpath('/home/ljp/Science/Loann/Loann Codes and Content/Scripts Matlab')


%% Define parameters, load cell and build info_fish

% parameters
fish = 3:18;
light = {' Off', ' On'};
rot = 0:2:12;

% cell
megapath = '/home/ljp/Science/Hippolyte/MegaFish.mat';
load(megapath);

% info_fish
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


%% Getting global regression using all rotations

% First isolate time series with at least two bout
keep_serie = (info_fish.info_time(:, 5) > 1);
time_keep = info_fish.info_time(keep_serie, :);
bout_keep = info_fish.info_bout(keep_serie, :);
keeplen = sum(keep_serie);

% - angle after bout:
outangle = {zeros(0, 1), zeros(0, 1)};
% - absolute angle before bout:
absangle = {zeros(0, 1), zeros(0, 1)};
% - absolute angle differential before bout:
absdiffang = {zeros(0, 1), zeros(0, 1)};
% - angle differential before bout:
diffang = {zeros(0, 1), zeros(0, 1)};
% - angle after bout on former bout
formout = {zeros(0, 1), zeros(0, 1)};
% NB: cells with left facing at index 1 and right facing at index 2


% Then get for each bout all the important information
for k = 1:keeplen
    % get right indexes
    time = time_keep(k, :);
    startpts = bout_keep{k, 1};
    endpts = bout_keep{k, 2};
    % get angle associated to try
    angle = MegaFish{time(1), time(2), time(3)}.angle;
    % loop over bouts
    for b = 2:length(endpts)
        if angle(startpts(b)) < 0 && -180 <= angle(startpts(b)) % -180 is because there is a weird point at -500 degrees
            left_or_right = 1; % left
        elseif 0 <= angle(startpts(b)) && angle(startpts(b)) <= 180
            left_or_right = 2; % right
        end
        outangle{left_or_right} = cat(1, outangle{left_or_right}, angle(endpts(b)));
        absangle{left_or_right} = cat(1, absangle{left_or_right}, abs(angle(startpts(b))));
        absdiffang{left_or_right} = cat(1, absdiffang{left_or_right}, abs(angle(startpts(b))-angle(startpts(b)-1)));
        diffang{left_or_right} = cat(1, diffang{left_or_right}, angle(startpts(b))-angle(startpts(b)-1));
        formout{left_or_right} = cat(1, formout{left_or_right}, angle(endpts(b-1)));
    end
end

figure
subplot(3, 2, 1)
hold on
plot(absangle{1}, outangle{1}, '.')
plot(absangle{2}, outangle{2}, '.')
legend('left facing', 'right facing')
title('angle after bout against absolute angle before bout', 'Interpreter', 'latex')
xlabel('angle [degrees]', 'Interpreter', 'latex')
ylabel('angle [degrees]', 'Interpreter', 'latex')
axis([-180, 180, -180, 180])
subplot(3, 2, 3)
hold on
plot(absdiffang{1}, outangle{1}, '.')
plot(absdiffang{2}, outangle{2}, '.')
legend('left facing', 'right facing')
title('angle after bout against absolute angle difference before bout', 'Interpreter', 'latex')
xlabel('angle [degrees]', 'Interpreter', 'latex')
ylabel('angle [degrees]', 'Interpreter', 'latex')
axis([0, 15, -180, 180])
subplot(3, 2, 5)
hold on
plot(diffang{1}, outangle{1}, '.')
plot(diffang{2}, outangle{2}, '.')
legend('left facing', 'right facing')
title('angle after bout against angle difference before bout', 'Interpreter', 'latex')
xlabel('angle [degrees]', 'Interpreter', 'latex')
ylabel('angle [degrees]', 'Interpreter', 'latex')
axis([-15, 15, -180, 180])
subplot(3, 2, 2)
hold on
plot(formout{1}, outangle{1}, '.')
plot(formout{2}, outangle{2}, '.')
legend('left facing', 'right facing')
title('angle after bout against angle after former bout', 'Interpreter', 'latex')
xlabel('angle [degrees]', 'Interpreter', 'latex')
ylabel('angle [degrees]', 'Interpreter', 'latex')
axis([-180, 180, -180, 180])
subplot(3, 2, 4)
hold on
plot(absangle{1}, formout{1}, '.')
plot(absangle{2}, formout{2}, '.')
legend('left facing', 'right facing')
title('angle after former bout against absolute angle before bout', 'Interpreter', 'latex')
xlabel('angle [degrees]', 'Interpreter', 'latex')
ylabel('angle [degrees]', 'Interpreter', 'latex')
axis([-180, 180, -180, 180])

% define test to train set proportion
test_prop = 0.1;

% define number of regressors and results cells
params = [1, 2, 2, 4, 5];
resids = cell(length(params), 2);
coeffs = cell(length(params), 2);
xtest = cell(length(params), 2);

% fill these cells
yout = cell(1, 2);
mix = cell(1, 2);
for lor = 1:2
    [coeffs{1, lor}, resids{1, lor}, ~, ~, xtest{1, lor}, yout{lor}, mix{lor}] = train_test_regress(test_prop, outangle{lor});
    [coeffs{2, lor}, resids{2, lor}, ~, ~, xtest{2, lor}] = train_test_regress(test_prop, outangle{lor}, absangle{lor}, mix{lor});
    [coeffs{3, lor}, resids{3, lor}, ~, ~, xtest{3, lor}] = train_test_regress(test_prop, outangle{lor}, formout{lor}, mix{lor});
    [coeffs{4, lor}, resids{4, lor}, ~, ~, xtest{4, lor}] = train_test_regress(test_prop, outangle{lor}, absangle{lor}, absdiffang{lor}, diffang{lor}, mix{lor});
    [coeffs{5, lor}, resids{5, lor}, ~, ~, xtest{5, lor}] = train_test_regress(test_prop, outangle{lor}, absangle{lor}, absdiffang{lor}, diffang{lor}, formout{lor}, mix{lor});
end

% compute matrix of F-tests
ftests = {zeros(length(params), length(params)), zeros(length(params), length(params))};
for lor = 1:2
    for i = 1:length(params)
        for j = 1:length(params)
            if j > i
                fstat = ((sum(resids{i, lor}.^2)-sum(resids{j, lor}.^2))/(j-i))/ (sum(resids{j, lor}.^2)/(length(yout{lor})-j));
                ftests{lor}(i, j) = fpdf(fstat, j-i, length(yout{lor})-j);
            end
        end
    end
end


%% Train regressions on a particular rotation and test on other

% First isolate time series with at least two bout
keep_serie = (info_fish.info_time(:, 5) > 1);
time_keep = info_fish.info_time(keep_serie, :);
bout_keep = info_fish.info_bout(keep_serie, :);
keeplen = sum(keep_serie);

% - angle after bout:
outangle = repmat({zeros(0, 1)}, length(rot), 2);
% - absolute angle before bout:
absangle = repmat({zeros(0, 1)}, length(rot), 2);
% - absolute angle differential before bout:
absdiffang = repmat({zeros(0, 1)}, length(rot), 2);
% - angle differential before bout:
diffang = repmat({zeros(0, 1)}, length(rot), 2);
% - angle after bout on former bout
formout = repmat({zeros(0, 1)}, length(rot), 2);
% NB: cells with left facing at index 1 and right facing at index 2

% save all information per rotation
for r = 1:length(rot)
    % get all time series associated to that rotation
    keep_rot = (time_keep(:, 3) == r);
    time_rot = time_keep(keep_rot, :);
    bout_rot = bout_keep(keep_rot, :);
    rotlen = sum(keep_rot);
    for k = 1:rotlen
        % get right indexes
        time = time_rot(k, :);
        startpts = bout_rot{k, 1};
        endpts = bout_rot{k, 2};
        % get angle associated to try
        angle = MegaFish{time(1), time(2), time(3)}.angle;
        % loop over bouts
        for b = 2:length(endpts)
            if angle(startpts(b)) < 0 && -180 <= angle(startpts(b)) % -180 is because there is a weird point at -500 degrees
                left_or_right = 1; % left
            elseif 0 <= angle(startpts(b)) && angle(startpts(b)) <= 180
                left_or_right = 2; % right
            end
            outangle{r, left_or_right} = cat(1, outangle{r, left_or_right}, angle(endpts(b)));
            absangle{r, left_or_right} = cat(1, absangle{r, left_or_right}, abs(angle(startpts(b))));
            absdiffang{r, left_or_right} = cat(1, absdiffang{r, left_or_right}, abs(angle(startpts(b))-angle(startpts(b)-1)));
            diffang{r, left_or_right} = cat(1, diffang{r, left_or_right}, angle(startpts(b))-angle(startpts(b)-1));
            formout{r, left_or_right} = cat(1, formout{r, left_or_right}, angle(endpts(b-1)));
        end
    end
end

figure
for r = 1:length(rot)
    subplot(2, length(rot), r)
    hold on
    plot(absangle{r, 1}, outangle{r, 1}, '.')
    plot(absangle{r, 2}, outangle{r, 2}, '.')
    legend('left facing', 'right facing')
    title(['x: angle before bout, rotation: ', num2str(rot(r))])
    xlabel('angle [degrees]', 'Interpreter', 'latex')
    ylabel('angle [degrees]', 'Interpreter', 'latex')
    axis([-180, 180, -180, 180])
    subplot(2, length(rot), length(rot)+r)
    hold on
    plot(formout{r, 1}, outangle{r, 1}, '.')
    plot(formout{r, 2}, outangle{r, 2}, '.')
    legend('left facing', 'right facing')
    title(['x: angle after former bout, rotation: ', num2str(rot(r))])
    xlabel('angle [degrees]', 'Interpreter', 'latex')
    ylabel('angle [degrees]', 'Interpreter', 'latex')
    axis([-180, 180, -180, 180])
end

% define test to train set proportion
test_prop = 0.2;

% get rmse for each rotation
rmse = {zeros(length(rot)), zeros(length(rot))};
for lor = 1:2
    for ri = 1:length(rot)
        % train test indexes
        boutlen = length(outangle{ri, lor});
        mix = randperm(boutlen)';
        test_num = round(test_prop * boutlen);
        train = mix(1:boutlen-test_num);
        test = mix(boutlen-test_num+1:end);
        % train multilinear regression
        [train_coeffs, ~, ~, mu_train, std_train] = train_regress(outangle{ri, lor}(train), absangle{ri, lor}(train), absdiffang{ri, lor}(train), ...
                                                                  diffang{ri, lor}(train), formout{ri, lor}(train));
        % get rmse on all rotations
        for rj = 1:length(rot)
            % if same rotation use test set
            if rj == ri
                test_residuals = test_regress(train_coeffs, outangle{rj, lor}(test), mu_train, std_train, absangle{rj, lor}(test), absdiffang{rj, lor}(test), ...
                                                                  diffang{rj, lor}(test), formout{rj, lor}(test));
            % else use all the data available
            else
                test_residuals = test_regress(train_coeffs, outangle{rj, lor}, mu_train, std_train, absangle{rj, lor}, absdiffang{rj, lor}, ...
                                                                  diffang{rj, lor}, formout{rj, lor});
            end
            % compute rmse from residuals
            rmse{lor}(ri, rj) = sqrt(mean(test_residuals.^2));
        end
    end
end

% average rmse and rmse transpose to have symmetry
for lor = 1:2
    rmse{lor} = (rmse{lor} + rmse{lor}')/2;
end

% plot
figure
left_or_right = {'left', 'right'};
for lor = 1:2
    subplot(3, 2, (lor-1)+[1, 3])
    image(rmse{lor}, 'CDataMapping', 'scaled')
    colorbar
    title(['RMSE for regression using one rotation on all rotations, fish facing ', left_or_right{lor}])
    xlabel('rotation [degrees]', 'Interpreter', 'latex')
    ylabel('rotation [degrees]', 'Interpreter', 'latex')
    xticklabels(rot)
    yticklabels(rot)
    subplot(3, 2, (lor-1)+5)
    plot(mean(rmse{lor}), 'd:')
    title(['mean RMSE per rotation, fish facing ', left_or_right{lor}])
    xlabel('rotation [degrees]', 'Interpreter', 'latex')
    ylabel('mean RMSE', 'Interpreter', 'latex')
end


%% Also separate light options

% First isolate time series with at least two bout
keep_serie = (info_fish.info_time(:, 5) > 1);
time_keep = info_fish.info_time(keep_serie, :);
bout_keep = info_fish.info_bout(keep_serie, :);
keeplen = sum(keep_serie);

% - angle after bout:
outangle = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - absolute angle before bout:
absangle = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - absolute angle differential before bout:
absdiffang = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - angle differential before bout:
diffang = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - angle after bout on former bout
formout = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% NB: cells with left facing at index 1 and right facing at index 2

% save all information per rotation
for r = 1:length(rot)
    for l = 1:length(light)
        % get all time series associated to that rotation and that light
        keep_rot = (time_keep(:, 3) == r & time_keep(:, 2) == l);
        time_rot = time_keep(keep_rot, :);
        bout_rot = bout_keep(keep_rot, :);
        rotlen = sum(keep_rot);
        for k = 1:rotlen
            % get right indexes
            time = time_rot(k, :);
            startpts = bout_rot{k, 1};
            endpts = bout_rot{k, 2};
            % get angle associated to try
            angle = MegaFish{time(1), time(2), time(3)}.angle;
            % loop over bouts
            for b = 2:length(endpts)
                if angle(startpts(b)) < 0 && -180 <= angle(startpts(b)) % -180 is because there is a weird point at -500 degrees
                    left_or_right = 1; % left
                elseif 0 <= angle(startpts(b)) && angle(startpts(b)) <= 180
                    left_or_right = 2; % right
                end
                outangle{r, l, left_or_right} = cat(1, outangle{r, l, left_or_right}, angle(endpts(b)));
                absangle{r, l, left_or_right} = cat(1, absangle{r, l, left_or_right}, abs(angle(startpts(b))));
                absdiffang{r, l, left_or_right} = cat(1, absdiffang{r, l, left_or_right}, abs(angle(startpts(b))-angle(startpts(b)-1)));
                diffang{r, l, left_or_right} = cat(1, diffang{r, l, left_or_right}, angle(startpts(b))-angle(startpts(b)-1));
                formout{r, l, left_or_right} = cat(1, formout{r, l, left_or_right}, angle(endpts(b-1)));
            end
        end
    end
end

% define test to train set proportion
test_prop = 0.2;

% get rmse for each rotation
rmse = {zeros(length(rot)*length(light)), zeros(length(rot)*length(light))};
for lor = 1:2
    for li = 1:length(light)
        for ri = 1:length(rot)
            % train test indexes
            boutlen = length(outangle{ri, li, lor});
            mix = randperm(boutlen)';
            test_num = round(test_prop * boutlen);
            train = mix(1:boutlen-test_num);
            test = mix(boutlen-test_num+1:end);
            % train multilinear regression
            [train_coeffs, ~, ~, mu_train, std_train] = train_regress(outangle{ri, li, lor}(train), absangle{ri, li, lor}(train), absdiffang{ri, li, lor}(train), ...
                                                                      diffang{ri, li, lor}(train), formout{ri, li, lor}(train));
            % get rmse on all rotations
            for lj = 1:length(light)
                for rj = 1:length(rot)
                    % if same rotation use test set
                    if rj == ri && lj == li
                        test_residuals = test_regress(train_coeffs, outangle{rj, lj, lor}(test), mu_train, std_train, absangle{rj, lj, lor}(test), absdiffang{rj, lj, lor}(test), ...
                                                                          diffang{rj, lj, lor}(test), formout{rj, lj, lor}(test));
                    % else use all the data available
                    else
                        test_residuals = test_regress(train_coeffs, outangle{rj, lj, lor}, mu_train, std_train, absangle{rj, lj, lor}, absdiffang{rj, lj, lor}, ...
                                                                          diffang{rj, lj, lor}, formout{rj, lj, lor});
                    end
                    % compute rmse from residuals
                    rmse{lor}((li-1)*length(rot)+ri, (lj-1)*length(rot)+rj) = sqrt(mean(test_residuals.^2));
                end
            end
        end
    end
end

% average rmse and rmse transpose to have symmetry
for lor = 1:2
    rmse{lor} = (rmse{lor} + rmse{lor}')/2;
end

% plot
figure
left_or_right = {'left', 'right'};
for lor = 1:2
    subplot(3, 2, (lor-1)+[1, 3])
    hold on
    image(rmse{lor}, 'CDataMapping', 'scaled')
    plot([7.5, 7.5], [0.5, 14.5], 'k')
    plot([0.5, 14.5], [7.5, 7.5], 'k')
    colorbar
    title(['RMSE for regression using one rotation on all rotations, fish facing ', left_or_right{lor}])
    xlabel('rotation [degrees]', 'Interpreter', 'latex')
    ylabel('rotation [degrees]', 'Interpreter', 'latex')
    rotnew = [rot; rot];
    xticklabels(rotnew(:))
    yticklabels(rotnew(:))
    subplot(3, 2, (lor-1)+5)
    hold on
    plot(mean(rmse{lor}), 'd:')
    plot([7.5, 7.5], [min(mean(rmse{lor})), max(mean(rmse{lor}))], 'k')
    title(['mean RMSE per rotation, fish facing ', left_or_right{lor}])
    xlabel('rotation [degrees]', 'Interpreter', 'latex')
    ylabel('mean RMSE', 'Interpreter', 'latex')
end


%% Scatter plots separating left and right

% First isolate time series with at least two bout
keep_serie = (info_fish.info_time(:, 5) > 1);
time_keep = info_fish.info_time(keep_serie, :);
bout_keep = info_fish.info_bout(keep_serie, :);
keeplen = sum(keep_serie);

% time increment
dt = MegaFish{1, 1, 1}.t(2) - MegaFish{1, 1, 1}.t(1);

% - tank angular speed
tankspeed = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - angle after bout:
corangle = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - absolute angle before bout:
absangle = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - absolute angle differential before bout:
absdiffang = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - angle differential before bout:
diffang = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - angle after bout on former bout
absformout = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - angular speed before bout fit on 5 points
fit5speed = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - angular speed before bout fit on 5 points
IBinterval = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% - angular speed before bout fit on 5 points
deltaang = repmat({zeros(0, 1)}, length(rot), length(light), 2);
% NB: cells with left facing at index 1 and right facing at index 2

% total matrices
X_fishtank = zeros(0, 8);
y_fishtank = zeros(0, 1);

% save all information per rotation
for r = 1:length(rot)
    for l = 1:length(light)
        % get all time series associated to that rotation and that light
        keep_rot = (time_keep(:, 3) == r & time_keep(:, 2) == l);
        time_rot = time_keep(keep_rot, :);
        bout_rot = bout_keep(keep_rot, :);
        rotlen = sum(keep_rot);
        for k = 1:rotlen
            % get right indexes
            time = time_rot(k, :);
            startpts = bout_rot{k, 1};
            endpts = bout_rot{k, 2};
            % get angle associated to try
            angle = MegaFish{time(1), time(2), time(3)}.angle;
            % loop over bouts
            for b = 2:length(endpts)
                if angle(startpts(b)) < 0 && -180 <= angle(startpts(b)) % -180 is because there is a weird point at -500 degrees
                    left_or_right = 1; % left
                elseif 0 <= angle(startpts(b)) && angle(startpts(b)) <= 180
                    left_or_right = 2; % right
                else
                    continue
                end
                % angle correction
                corangle{r, l, left_or_right} = cat(1, corangle{r, l, left_or_right}, abs(angle(startpts(b)))-abs(angle(endpts(b))));
                % absolute angle before bout
                absangle{r, l, left_or_right} = cat(1, absangle{r, l, left_or_right}, -abs(angle(startpts(b))));
                % absolute angle difference and angle difference before bout
                absdiffang{r, l, left_or_right} = cat(1, absdiffang{r, l, left_or_right}, abs(angle(startpts(b))-angle(startpts(b)-1))/dt);
                diffang{r, l, left_or_right} = cat(1, diffang{r, l, left_or_right}, (abs(angle(startpts(b)))-abs(angle(startpts(b)-1)))/dt);
                % absolute angle after former bout
                absformout{r, l, left_or_right} = cat(1, absformout{r, l, left_or_right}, abs(angle(endpts(b-1))));
                % fit on 10 points before bout to have instantaneous speed
                pts_5 = -abs(angle(startpts(b)-4:startpts(b)))';
                xtrain = [ones(5, 1), (1:5)'];
                coeffs = (xtrain' * xtrain) \ xtrain' * pts_5;
                fit5speed{r, l, left_or_right} = cat(1, fit5speed{r, l, left_or_right}, coeffs(2)/dt);
                % interbout interval and angular speed from last bout
                IBI = (startpts(b) - endpts(b-1)) * dt;
                IBinterval{r, l, left_or_right} = cat(1, IBinterval{r, l, left_or_right}, IBI);
                deltaang{r, l, left_or_right} = cat(1, deltaang{r, l, left_or_right}, (abs(angle(endpts(b-1)))-abs(angle(startpts(b))))/IBI);
                % tank speed
                tankspeed{r, l, left_or_right} = cat(1, tankspeed{r, l, left_or_right}, rot(r));
            end
        end
    end
end


%% Scatter plots

% legends
lor_legend = {'left facing', 'right facing'};
r_legend = {'0 degrees', '2 degrees', '4 degrees', '6 degrees', '8 degrees', '10 degrees', '12 degrees'};
l_legend = {'lights off', 'lights on'};

% what to plot
xplot = absformout;
yplot = corangle;

% emphasis on fish facing
figure
hold on
for lor = 1:2
    xtemp = zeros(0, 1);
    ytemp = zeros(0, 1);
    for r = 1:length(rot)
        for l = 1:length(light)
            xtemp = cat(1, xtemp, xplot{r, l, lor});
            ytemp = cat(1, ytemp, yplot{r, l, lor});
        end
    end
    scatter(xtemp, ytemp, 75, '.')
end
legend(lor_legend)
title('correction angle against absolute angle before bout')
axis equal

% emphasis on rotation
figure
hold on
for r = 1:length(rot)
    xtemp = zeros(0, 1);
    ytemp = zeros(0, 1);
    for lor = 1:2
        for l = 1:length(light)
            xtemp = cat(1, xtemp, xplot{r, l, lor});
            ytemp = cat(1, ytemp, yplot{r, l, lor});
        end
    end
    scatter(xtemp, ytemp, 75, '.')
end
legend(r_legend)
title('correction angle against absolute angle before bout')
axis equal

% emphasis on light
figure
hold on
for l = 1:length(light)
    xtemp = zeros(0, 1);
    ytemp = zeros(0, 1);
    for lor = 1:2
        for r = 1:length(rot)
            xtemp = cat(1, xtemp, xplot{r, l, lor});
            ytemp = cat(1, ytemp, yplot{r, l, lor});
        end
    end
    scatter(xtemp, ytemp, 75, '.')
end
legend(l_legend)
title('correction angle against absolute angle before bout')
% axis equal


%% Saving all important data in matrices

% First isolate time series with at least two bout
keep_serie = (info_fish.info_time(:, 5) > 1);
time_keep = info_fish.info_time(keep_serie, :);
bout_keep = info_fish.info_bout(keep_serie, :);
keeplen = sum(keep_serie);

% time increment
dt = MegaFish{1, 1, 1}.t(2) - MegaFish{1, 1, 1}.t(1);

% total matrices
X_fishtank = zeros(0, 8);
y_fishtank = zeros(0, 1);

% Then get for each bout all the important information
for k = 1:keeplen
    % get right indexes
    time = time_keep(k, :);
    startpts = bout_keep{k, 1};
    endpts = bout_keep{k, 2};
    % get angle associated to try
    angle = MegaFish{time(1), time(2), time(3)}.angle;
    % loop over bouts
    for b = 2:length(endpts)
        if angle(startpts(b)) < 0 && -180 <= angle(startpts(b)) % -180 is because there is a weird point at -500 degrees
            lor = -1; % left
        elseif 0 <= angle(startpts(b)) && angle(startpts(b)) <= 180
            lor = 1; % right
        else
            continue
        end        
        % input X_fishtank
        lightstatus = (time(2)-1.5) * 2; % light status
        tankspeed = rot(time(3)); % tank rotation speed
        absangle = -abs(angle(startpts(b))); % minus absolute angle before bout
        absformout = abs(angle(endpts(b-1))); % absolute angle after former bout
        pts_5 = -abs(angle(startpts(b)-4:startpts(b)))';
        xtrain = [ones(5, 1), (1:5)'];
        coeffs = (xtrain' * xtrain) \ xtrain' * pts_5;
        fit5speed = coeffs(2)/dt; % fit on last 5 points for speed
        IBI = (startpts(b) - endpts(b-1)) * dt; % interbout interval
        deltaang = (abs(angle(endpts(b-1)))-abs(angle(startpts(b))))/IBI; % angle difference from last bout divided by IBI
        % total matrix X_fishtank
        X_fishtank = cat(1, X_fishtank, [lightstatus, tankspeed, lor, absangle, absformout, fit5speed, IBI, deltaang]);
        % output y_fishtank
        y_fishtank = cat(1, y_fishtank, abs(angle(startpts(b)))-abs(angle(endpts(b))));
    end
end








