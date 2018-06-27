function [hValues, pValues , comp] = statDiffAnalysis(params,pMax)
% [hValues, pValues , comp] = statDiffAnalysis(params,pMax)
%
% This function performs a series of Mann-Whitney U-tests to determine if
% the parameters passed to the function show statistically significant
% difference between them.
%
% 'params' is a cell array of size NxM, where N is the number of channels
% and M is the number time points of the signals being analysed. Every
% entry is the experimentally calibrated population of models for channel n
% at timestep m. All ePOMs are tested for statistically significant
% difference with the ePOM provided in params{:,1}.
%
% 'pMax' is the level of statistical significance to be used in the
% analysis. It must be a number between 0 and 1 representing the maximum
% accepted probability that two means come from the same sample.
%
% 'comp' is a cell array of 1xlength(pMax). Each entry shows which
% comparisons were made to find the values shown in 'hValues' and 'pValues'
%
% 'hValues' is a cell array of 1xlength(pMax). Each entry is a cell array
% containing the results of doing the comparisons shown in 'comp' on every
% channel of params. On each entry will contain a 1 if the null hypothesis
% was rejected for a given parameter and 0 otherwise.
%
% 'pValues' has the same structure as hValues, but returns the probability
% that the means compared come from the same sample.
%
% This function was developed in the Multiscale Cardiovascular Engineering
% Group (MUSE) at University College London by Carlos Ledezma.

% Number of channels in the signals
numChans = size(params,1);
% Number of time points
numSigs = size(params,2);
% Number of parameters in the model
numParams = size(params{1},2);
% Number of t-tests that need to be performed (combination of numSigs taken
% in pairs)
numTtest = length(numSigs)-1; 

% Initialise output variables
pValues = cell(numChans,numTtest);
hValues = pValues;
comp = cell(1,numTtest);
for j = 2 : numSigs
    comp{j-1} = ['1-' num2str(j)];
end

for i = 1 : numChans % For every channel
    t = 1;
    % Comparison is made sequentially (1-1, 1-2,..., 1-N)
    thisSig = params{i,1};% Current PoM
    if ~isnan(params{i,1}(1))
        for k = 2 : numSigs
            nextSig = params{i,k};% PoM to be compared with
            if ~isnan(params{i,k}(1))
                for l = 1 : numParams % Perform MWW U-test for every parameter
                    %Mann-Whitney U-test
                    [pValues{i,t}(l), hValues{i,t}(l)] = ranksum(thisSig(:,l), nextSig(:,l), 'alpha', pMax);
                end
            end
            t = t + 1;
        end
    end
end