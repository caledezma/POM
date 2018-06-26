function [params, ind] = findParamsFromPoM(PoM,ARI,dt,CL)

% Initialise output parameters
params = cell(length(ARI),1);
ind = params;

CL = find(PoM(1).CL == CL);

for i = 1 : length(ARI) % Iterate over every channel
    thisARI = ARI{i}(:,3)*dt; %Recover the ARI of the current channels
    % Define the range of tolerance for the calibration
    minAPD = (mean(thisARI)-std(thisARI))*1000;
    maxAPD = (mean(thisARI)+std(thisARI))*1000;
    
    k = 1;
    % Calibrate the population
    for j = 1 : length(PoM) % Iterate over every model
        % If the ARI of a given model falls within the accepted range, save
        % the parameters and the index.
        if PoM(j).APD(CL) >= minAPD && PoM(j).APD(CL) <= maxAPD
            params{i,1}(k,:) = PoM(j).params;
            ind{i,1}(k) = j;
            k=k+1;
        end
    end
    
    % If no model in the population falls within the acceptance range,
    % return NaN.
    if k == 1 
        params{i,1} = NaN*ones(1,length(PoM(1).params));
    end
end