clear; close all; clc
addpath('/home/ljp/Science/Loann/Loann Codes and Content/Scripts Matlab')


%% Load files of parameters

% Load inputs
X_ft = readtable('X_fishtank.csv', 'ReadVariableNames', false);
X_ft = table2array(X_ft);

% load outputs
y_ft = readtable('y_fishtank.csv', 'ReadVariableNames', false);
y_ft = table2array(y_ft);

% reminder of the input labels
labels = {'light', 'rotation', 'orientation', 'absolute angle before bout', ...
          'absolute angle after last bout', 'speed fitted on 5 last points', ...
          'interbout interval', 'angle difference from last bout / IBI'};


%% Shuffle and split into training and test set

% important parameters
[m, n] = size(X_ft);

% train and test vectors
test_ratio = 0.1;
shuffle = randperm(m);
bval = round((1-test_ratio)*m);
train = shuffle(1:bval);
test = shuffle(bval+1:end);

% get training and test sets
xtrain = X_ft(train, :);
ytrain = y_ft(train, :);
xtest = X_ft(test, :);
ytest = y_ft(test);


%% Constant model

% regressors training
regtrain = ones(bval, 1);
% regressors test
regtest = ones(m-bval, 1);

% coefficients
coeffs = (regtrain' * regtrain) \ regtrain' * ytrain;

% residuals
resids_null = ytest - regtest * coeffs;

% RMS of residuals
RMSr = sqrt(mean(resids_null.^2));

% print information
fprintf('Root mean square of residuals for constant model is %.2f \n', RMSr);


%% Schoppik model
% not exactly Schoppik model but, the most similar possible with my params

% parameters to keep
park = [xtrain(:, [3, 4, 6]), abs(xtrain(:, 6))];
parktest = [xtest(:, [3, 4, 6]), abs(xtest(:, 6))];

% regressors training
nr = size(park, 2);
mu_train = zeros(1, nr);
std_train = zeros(1, nr);
regtrain = ones(bval, nr+1);
for i = 1:nr
    mu_train(i) = mean(park(:, i));
    std_train(i) = std(park(:, i));
    regtrain(:, i+1) = (park(:, i) - mu_train(i)) ./ std_train(i);
end
% regressors test
regtest = ones(m-bval, nr+1);
for i = 1:nr
    regtest(:, i+1) = (parktest(:, i) - mu_train(i)) ./ std_train(i);
end

% coefficients
coeffs = (regtrain' * regtrain) \ regtrain' * ytrain;

% residuals
resids = ytest - regtest * coeffs;

% RMS of residuals
RMSr = sqrt(mean(resids.^2));

% F-statistic
RSS1 = sum(resids_null.^2);
RSS2 = sum(resids.^2);
fstat = ((RSS1-RSS2)/((nr+1)-1)) / (RSS2/((m-bval)-(nr+1)));
ftest = fpdf(fstat, (nr+1)-1, (m-bval)-(nr+1));

% print information
fprintf('Root mean square of residuals for constant model is %.2f \n', RMSr);
fprintf('Associated F-statistic is %.2f, ftest associated is %.2f \n', [fstat, ftest]);


%% Full model
% model using all the parameters available

% parameters to keep
park = xtrain;
parktest = xtest;

% regressors training
nr = size(park, 2);
mu_train = zeros(1, nr);
std_train = zeros(1, nr);
regtrain = ones(bval, nr+1);
for i = 1:nr
    mu_train(i) = mean(park(:, i));
    std_train(i) = std(park(:, i));
    regtrain(:, i+1) = (park(:, i) - mu_train(i)) ./ std_train(i);
end
% regressors test
regtest = ones(m-bval, nr+1);
for i = 1:nr
    regtest(:, i+1) = (parktest(:, i) - mu_train(i)) ./ std_train(i);
end

% coefficients
coeffs = (regtrain' * regtrain) \ regtrain' * ytrain;

% residuals
resids = ytest - regtest * coeffs;

% RMS of residuals
RMSr = sqrt(mean(resids.^2));

% F-statistic
RSS1 = sum(resids_null.^2);
RSS2 = sum(resids.^2);
fstat = ((RSS1-RSS2)/((nr+1)-1)) / (RSS2/((m-bval)-(nr+1)));
ftest = fpdf(fstat, (nr+1)-1, (m-bval)-(nr+1));

% print information
fprintf('Root mean square of residuals for constant model is %.2f \n', RMSr);
fprintf('Associated F-statistic is %.2f, ftest associated is %.2f \n', [fstat, ftest]);










