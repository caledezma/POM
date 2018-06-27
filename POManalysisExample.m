% This script is an example of how to create a POM, calibrate it to some
% experimental data and obtain results from the calibration process. The
% script uses the TP06 model to create the POM, only a small population
% (1000 models) is created at a single cycle length (600 ms). This POM is
% then calibrated to three different experimental data:
% 'ARIdata/time_point_1.mat', 'ARIdata/time_point_2.mat',
% 'ARIdata/time_point_3.mat'. The population is saved in 'POM.mat' and the
% ePOMs are saved in 'ePOMs/ePOM1.mat', 'ePOMs/ePOM2.mat',
% 'ePOMs/ePOM3.mat'.
%
% Several outputs are saved in 'Results/':
%   'Boxplots/chanX.png' contains the boxplots of the statistical
%   distributions of the ePOMs' parameters. Each image contains the results
%   of calibrating the POM to channel X at of the three time points.
%
% This script was developed in the Multiscale Cardiovascular Engineering 
% Group at University College London (MUSE-UCL) by Carlos Ledezma. It's 
% protected by a Creative Commons Attribution-ShareAlike 4.0 International 
% license (https://creativecommons.org/licenses/by-sa/4.0/).
%
% For further details please refer to the published work. This paper must
% be cited if using any of the routines or data provided with this script:
% C. Ledezma et al. "Bridging organ and cellular-level behavior in ex-vivo
% experimental platforms using populations of models of cardiac
% electrophysiology" (2018). ASME Journal of Engineering and Science in
% Medical Diagnostics and Therapy. doi: 10.1115/1.4040589.

addpath('func');

% Names of TP06 parameters
paramNames = {'Gna','Gto','Gks','Gkr','Gk1','Pnak','Gcal','knaca','Gpca','Gpk','Gbna','Gbca'};

% Significance level to use in the statistical analysis
pMax = 0.001;

% Name for the population of models
POMpath = 'POM.mat';
% Paths to the ARI signals
ARIpath = {'ARIdata/time_point_1.mat',....
           'ARIdata/time_point_2.mat',...
           'ARIdata/time_point_3.mat' };
       
% Create folder to save the ePOMs
ePOMroot = 'ePOMs';
if isempty(dir(ePOMroot))
    mkdir(ePOMroot);
end
% Paths to the ePOMs
ePOMpath = {'ePOMs/ePOM1.mat', 'ePOMs/ePOM2.mat', 'ePOMs/ePOM3.mat'};

% Create folders to save the results
if isempty(dir('Results'))
    mkdir('Results')
end
if isempty(dir('Results/Boxplots'))
    mkdir('Results/Boxplots')
end
if isempty(dir(['Results/MWUtestResults' num2str(pMax*100)]))
    mkdir(['Results/MWUtestResults' num2str(pMax*100)]);
end

% Create a population of TP06 models and save it
% Note that this is just an example, for the parameters required to have a
% sufficiently larged, converged, population refer to the paper
POM = generatePoM(2,600,30,1000,0,2);
save(POMpath,'POM');

% Calibrate pooulation to each of the ARI data
ePOM1 = calibratePoM(POMpath,ARIpath{1});
save(ePOMpath{1},'ePOM1')
ePOM2 = calibratePoM(POMpath,ARIpath{2});
save(ePOMpath{2},'ePOM2')
ePOM3 = calibratePoM(POMpath,ARIpath{3});
save(ePOMpath{3},'ePOM3')

% Save the boxplots of the ePOMs' parameter distributions per channel
chanNum = length(ePOM1.params); % Number of channels in data
% Create the figure so that it runs in the background
monPos = get(0,'MonitorPositions');
fig = figure('outerposition',monPos(2,:));
set(fig,'visible','off');
set(fig,'PaperPositionMode','auto');
% Initialise the waitbar
wb = waitbar(0,'Generating boxplots');
% Labels for the boxplots
labels{1} = ARIpath{1}(9:end-4);
labels{2} = ARIpath{2}(9:end-4);
labels{3} = ARIpath{3}(9:end-4);

for i = 1 : chanNum % Loop over all the channels
    % Put all parameters in a single variable
    allParams = [ePOM1.params{i} ; ePOM2.params{i} ; ePOM3.params{i}];
    % Create groups for boxplots
    G = [ones(size(ePOM1.params{i},1),1) ; ...
        2*ones(size(ePOM2.params{i},1),1) ;...
        3*ones(size(ePOM3.params{i},1),1) ];
    for k = 1 : size(allParams,2) % Boxplot each parameter in a different plot
        subplot(3,4,k)
        boxplot(allParams(:,k), G, 'Labels', labels);
        title(['Channel ' num2str(i) , ', Parameter: ' paramNames{k}]);
    end
    waitbar(i/chanNum,wb,['Generating boxplots... ' num2str(i) '/' num2str(chanNum)]);
    saveas(fig,['Results/Boxplots/chan' num2str(i) '.png']);
end
close(fig);
close(wb);

%Temporarily turn off the warning of image too large for page
warnId = 'MATLAB:print:FigureTooLargeForPage';
warning('off',warnId);

% Perform statistical analysis (Mann-Whitney-Wilcoxon U-Test) comparing
% ePOM2 and ePOM3 to ePOM1 in search for significant differences

[hValues, pValues, comp] = ...
    statDiffAnalysis([ePOM1.params,ePOM2.params,ePOM3.params],pMax);

%Plot the significant changes and save them in the set directory
for i = 1 : size(hValues,2)
    % chanCount counts how many channels showed statistical
    % significance for each of the parameters.
    chanCount{i} = zeros(1,length(paramNames));
    % Set figure to run in the background
    fig = figure('outerposition',monPos(2,:));
    set(fig,'visible','off');
    set(fig,'PaperPositionMode','auto');
    hold on;
    % Create plot with marks on the channels that showed statistical
    % significance
    for j = 1 : size(hValues,1)
        for l = 1 : length(hValues{j,i})
            if hValues{j,i}(l)
                mark = '*b';
                plot(l,j,mark)
                chanCount{i}(l) = chanCount{i}(l) + 1;
            end
        end
    end
    % Ploting configurations
    ax = gca;
    ax.XTick = 1:length(hValues{1,1});
    ax.XTickLabel = paramNames;
    set(ax,'fontsize',18)
    xlim([0 length(hValues{1,1})+1]);
    ylim([0 size(hValues,1)+1])
    ylabel Channel;
    title(['When comparing ' comp{i}]);
    legend(['\rho \leq ' num2str(pMax)]);
    saveas(fig,['Results/MWUtestResults' num2str(pMax*100) '/comp' comp{i} '.png']);
    close(fig);
end

% Turn on warning again
warning('on',warnId);
