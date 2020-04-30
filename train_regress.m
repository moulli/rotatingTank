%% Function for multilinear regression

function [train_coeffs, X_train, y_train, mu_train, std_train] = train_regress(output, varargin)
% returns coefficients from training set.
    
    % size of output
    boutlen = length(output);
    
    % all outputs
    y_train = reshape(output, boutlen, 1);
    
    % define mu_train and std_train
    lenv = length(varargin);
    mu_train = zeros(lenv, 1);
    std_train = zeros(lenv, 1);
    
    % all regressors
    X_train = ones(boutlen, length(varargin)+1);
    for i = 1:lenv
        vari = reshape(varargin{i}, boutlen, 1);
        mu_train(i) = mean(vari);
        std_train(i) = std(vari);
        X_train(:, i+1) = (vari - mean(vari)) ./ std(vari);
    end
    
    % compute regressors
    train_coeffs = (X_train' * X_train) \ X_train' * y_train;

end