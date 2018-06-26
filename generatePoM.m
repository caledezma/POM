function PoM = generatePoM(model, CL, numBeats, popSize, lowMult, highMult)
% PoM = generatePoM(model, CL, numBeats, popSize)
%
% Create a population of models using latin hypercube sampling. 
%
% Inputs:
%   model = 1 (ORd) or 2 (TP06)
%   CL: array containing the cycle lengths at which each model in the
%       population will be stimulated
%   numBeats: number of beats needed to reach steady state
%   popSize: number of models in the population
%   lowMult: lower multiplier for the mean value of the parameters
%   highMult: higher multiplier for the mean value of the parameters
%
% Output:
%   PoM: a struct array containing all the models from the population that
%   were physiological. The struct contains the following fields:
%       V: a cell array with the action potential traces at each cycle
%       lenght
%       t: a cell array with the time vectors at each cycle lenght
%       params: array containing the values of the parameters used to solve 
%       the model that produced the AP trace
%       APD: the APD90 measured from each action potential trace
%       CL: the cycle lengths at which the model was solved
%
% This function was developed in the Multiscale Cardiovascular Engineering 
% Group at University College London (MUSE-UCL) by Carlos Ledezma. It's 
% protected by a Creative Commons Attribution-ShareAlike 4.0 International 
% license (https://creativecommons.org/licenses/by-sa/4.0/)


addpath('models','func')

% Choose variables depending on the model to be solved
if model == 1
    mod = @myORd;
    cellType = 1;
    meanParams = conductancesORd(cellType);
    savePath = 'myPoM/';
elseif model == 2 
    mod = @myTP06;
    cellType = 3;
    meanParams = conductancesTP06(cellType);
    savePath = 'myPoM/';
else
    error('Model has to be 1 (ORd) or 2 (TP06)');
end

% Sample the parameter space to obtain the PoM
params = generateSamples(popSize, meanParams, lowMult*meanParams, highMult*meanParams ,...
    popSize, 2, []);

numIt = size(params,1);

for i = 1 : numIt
    result(i) = struct('V',{{0}},'time',{{0}},'APD',[]);
end

if isempty(gcp('nocreate'))
    parpool();
end

% Temporary directory to save PoM pieces
if isempty(dir('myPoM'))
    mkdir('myPoM')
end

parfor_progress(numIt);
% Solve model for every combination of parameters
parfor i = 1 : numIt
    for j = 1 : length(CL)
        if model == 1
            X0 = initCondsORd;
        elseif model == 2
            X0 = initCondsTP06();
        else
            error('Model has to be 1 (ORd) or 2 (TP06)');
        end
        
        options = [];%odeset('MaxStep',1,'RelTol',1e-7,'AbsTol',1e-9);
        for n=1:numBeats
            [result(i).time{j} , X]=ode15s(mod,[0 CL(j)],X0,options,1,params(i,:),cellType);
            X0=X(size(X,1),:);
        end
        
        result(i).V{j} = X(:,1);
        result(i).APD(j) = findAPD(result(i).V{j},result(i).time{j},90);
               
    end
    
    % Save each model
    V = result(i).V;
    APD = result(i).APD;
    time = result(i).time;
    parsave([savePath 'PoM' num2str(i) '.mat'],{V , time , params(i,:) , APD, CL});
    parfor_progress;
end
parfor_progress(0);

clear result

disp('Saving population to a single variable...')
PoM = struct('V',[],'t',[],'params',[],'APD',[],'CL',[]);

% After creating the population, join all the pieces. Only include those
% that are physiological

k = 1;
for i = 1 : numIt
    allPhysio = true;
    load(['myPoM/PoM' num2str(i) '.mat']); % Load the record
    for j = 1 : length(CL) % Loop through the cycle lengths to check that all are physiological
        if ~isPhysio(var{1}{j},var{2}{j},var{4}(j))
            allPhysio = false;
        end
    end
    
    if allPhysio % If everything is in order, save the action potentials
        PoM(k).V = var{1};
        PoM(k).t = var{2};
        PoM(k).params = var{3};
        PoM(k).APD = var{4};
        PoM(k).CL = var{5};
        k = k + 1;
    end
    delete(['myPoM/PoM' num2str(i) '.mat']);
end

rmdir('myPoM')
