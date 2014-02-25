function cost = comparemodels(expfeatures, simfeatures, searchspace)

%{
PSEUDOCODE:
% get feature means and standard deviations from expfeatures
% ...
% compute experimentally-normalized simulated feature z-scores (for each simulation; i.e., N per model/simfile)
% ...
% compute mean simulated feature z-scores (one per model <=> one per point in searchspace)
% ...
% determine # params for each model
% ...
% compute model cost = f(feature z-scores, # params)
% ...
%}






%% COPIED FROM DEV WORK IN CharacterizeCells()

feature_weights = ones(size(features));

% calculate experimental mean and standard deviations per feature
% ...

% for tests set to 0s and 1s.
expt_mu = zeros(size(features));
expt_sd = ones(size(features));

% calculate z-scores on simulated feature vectors
zfeatures = (features - expt_mu) ./ expt_sd;

% calculate mean square feature z-scores (MSFZ)
MSFZ = mean(zfeatures.^2);
MSWFZ = mean((zfeatures.*feature_weights).^2); % weighted version

