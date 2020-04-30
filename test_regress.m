%% Function for multilinear regression

function [test_residuals, X_test, y_test] = test_regress(train_coeffs, output, mu_train, std_train, varargin)
% returns coefficients from training set and residuals from test set.
% training set and test sets are computed randomly using test_prop: the
% proportion of examples in the test set. varargin is a cell of regressors.
% if the last element of varargin is a vector of integer ranging from 1 to
% length(output), then it will be used as the mix vector
    
    % size of output
    boutlen = length(output);
    
    % all outputs
    y_test = reshape(output, boutlen, 1);
    
    % all regressors
    X_test = ones(boutlen, length(varargin)+1);
    for i = 1:length(varargin)
        X_test(:, i+1) = (varargin{i} - mu_train(i)) ./ std_train(i);
    end
    
    % compute residuals
    test_residuals = y_test - X_test * train_coeffs;

end