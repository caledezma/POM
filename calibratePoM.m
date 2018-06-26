function ePOM = calibratePoM(POMfile, dataFile)
%
% ePOM = calibratePoM(POMfile, dataFile)
% 
% Calibrate a previously generated Population of Models (POM) to
% experimental data to create an experimentally-calibrated population of
% models(ePOM)
%
% Inputs:
%   POMfile: a string containing the path to a .mat file containing a
%   previously constructed population of models. Preferably, this
%   population should be created using the 'generatePoM.m' routine provided
%   with this repository. The variable stored in 'POMfile' MUST be called
%   POM.
%   dataFile: a string containing the path to the activation-recovery
%   intervals (ARI) to which the population of models will be calibrated.
%   This datafile MUST contain two variables:
%       fs: a double specifying the sampling frequency at which the data
%       was acquired.
%       CL: a double specifying the cycle length at wich the data was
%       acquired.
%       ARI: a cell array where each entry ARI{i} corresponds to data
%       acquired by a single channel. ARI{i} must contain a double matrix,
%       each row of this matrix must be of the form:
%       ARI{i}(j,:) = [ActivationTime, RecoveryTime, ARI]
%
% Output:
%   ePOM: population containing only the models from the POM that
%   reproduced the behaviour observed in the data provided. ePOM is a
%   struct containing two fields:
%       params: is a cell array of the same size as ARI. Each entry
%       params{i} is a matrix where each row contains the paramters from a
%       model of the POM that reproduces the behaviour observed in ARI{i}.
%       ind: is a cell array of the same size as ARI. Each entry ind{i}
%       is a vector of length size(params{i},1), where ind(j) is the index
%       in the original POM where the model of parameters param{i}(j,:) can
%       be found.
%
% This function was developed in the Multiscale Cardiovascular Engineering 
% Group at University College London (MUSE-UCL) by Carlos Ledezma. It's 
% protected by a Creative Commons Attribution-ShareAlike 4.0 International 
% license (https://creativecommons.org/licenses/by-sa/4.0/)

addpath('func')

% Load the pre-calculated population of models
load(POMfile);

% The AP and time variables are not needed and take too much memory
POM = rmfield(POM,'V');
POM = rmfield(POM,'t');

% Load the ARI signals and sampling frequency
load(dataFile)

% Calibrate the POM to each of the provided channels
[ePOM.params, ePOM.ind] = findParamsFromPoM(POM, ARI, 1/fs, CL);


