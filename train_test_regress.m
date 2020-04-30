%% Function for multilinear regression

function [train_coeffs, test_residuals, X_train, y_train, X_test, y_test, mix] = train_test_regress(test_prop, output, varargin)
% returns coefficients from training set and residuals from test set.
% training set and test sets are computed randomly using test_prop: the
% proportion of examples in the test set. varargin is a cell of regressors.
% if the last element of varargin is a vector of integer ranging from 1 to
% length(output), then it will be used as the mix vector
    
    % size of output
    boutlen = length(output);

    % create training and test indexes
    if ~isempty(varargin)
        potential_mix = varargin{end};
        potential_mix = reshape(potential_mix, length(potential_mix), 1);
        if isequal(sort(potential_mix), (1:boutlen)')
            mix = potential_mix;
            varargin = varargin(1:end-1);
        else
            mix = randperm(boutlen)';
        end
    else
        mix = randperm(boutlen)';
    end
    test_num = round(test_prop * boutlen);
    train = mix(1:boutlen-test_num);
    test = mix(boutlen-test_num+1:end);
    
    % all outputs
    y_all = reshape(output, boutlen, 1);
    y_train = y_all(train);
    y_test = y_all(test);
    
    % all regressors
    X_all = ones(boutlen, length(varargin)+1);
    for i = 1:length(varargin)
        X_all(:, i+1) = (varargin{i} - mean(varargin{i})) ./ std(varargin{i});
    end
    X_train = X_all(train, :);
    X_test = X_all(test, :);
    
    % compute regressors
    train_coeffs = (X_train' * X_train) \ X_train' * y_train;
    
    % compute residuals
    test_residuals = y_test - X_test * train_coeffs;

end