function samples = generateSamples(numSamples, means, variL , variH , popSize, sampMethod, prevParams)

samples = zeros(numSamples,length(means));
population = zeros(popSize, length(means));
numParams = length(means);

%Build the population of parameters
for i = 1 : numParams
    step = (variH(i)-variL(i))/popSize;
    population(:,i) = variL(i) : step : variH(i) - step;
end

%Perform the sampling of the population
switch sampMethod
    case 1 %Uniform distribution sampling
        for i = 1 : numParams
            %Generate random samples following a uniform distribution in
            %the interval [1, popSize]
            randSamples = floor(1 + (popSize-1)*rand(numSamples,1));
            %Save those samples
            samples(:,i) = population(randSamples,i);
        end
        
    case 2 %Latin hypercube sampling
        if popSize ~= numSamples
            error('For proper LHS the size of the population must be the same as the number of samples to be taken');
        end
        
        if ~isempty(prevParams)
            % Eliminate previously taken parameters
            for i = 1 : size(prevParams,1)
                for j = 1 : size(prevParams,2)
                    population(population(:,j) == prevParams(i,j),j) = -1;
                end
            end
            
            popSize = popSize-size(prevParams,1);
            population(population == -1) = [];
            population = reshape(population,popSize,numParams);
        end
        
        % Reshape the solution arrays
        numSamples = popSize;
        samples = zeros(numSamples,length(means));
        
        for i = 1 : numSamples
            for j = 1 : numParams
                %Pick a random value for a given parameter (uniform
                %distribution)
                randSample = floor(1 + (popSize-1)*rand(1));
                %Save the chosen parameter
                samples(i,j) = population(randSample,j);
                %Mark parameter as taken
                population(randSample,j) = -1;
            end
            %Eliminate taken samples from population
            population(population==-1) = [];
            popSize = popSize-1;
            population = reshape(population,popSize,numParams);
        end
        
    case 3 %Box and Muller sampling
        for i = 1 : numParams
            %Generate the random points
            randSamples = sqrt(-2*log(rand(numSamples,1))).*cos(2*pi*rand(numSamples,1));
            randSamples = round(linTrans(randSamples,popSize,1));
            samples(:,i) = population(randSamples,i);
        end
end