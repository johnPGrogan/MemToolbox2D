% Demo2D_1.m
clear all
% demo code from the Tutorial.pdf


%% fitting to real data

% your data structure needs at minimum the data.errors field, which is
% response location minus target location.
% If you model has misbinding, you also need data.distractors, which is the
% distractors locations minus the target locations.
% 
% A default screen size of [1366; 768] will be assumed unless you specify
% the data.dimensions. You should change this to match your screen
% 
% Other models may require extra fields, such as data.resps or data.targets
% or data.distractors1 and data.distractors2.
% 
% 

%% set up trials

numTrials = 300; % number of trials
itemsPerTrial = ones(1,numTrials)*3; % 3 items per trial
dimensions = [1366; 768]; % x and y dimensions in units (pixels here)
simulatedData = GenerateDisplays2D(numTrials, itemsPerTrial, 1,dimensions); % make locations, pick target, calculate target-distractor distance
% distractors = distractorLocation - targetLocation

% simulate responses

model = SwapModel2D(); % the model to simulate
params = [0.1, 0.1, 20]; % gamma, beta, SD

simulatedData.errors = SampleFromModel2D(model, params,...
                 [1 numTrials], simulatedData); 
% simulate model to generate responses, calculate distance to targets
% errors = responseLocation - targetLocation

%% fit the data

fit = MLE(simulatedData, model); % fit the model to simulated data
% maximum likelihood method
% parameters are in order of model.paramNames

% you can also use Bayesian MCMC fitting (takes longer)
fit2 = MemFit2D(simulatedData, model);
% 'Verbosity',0 prevents it asking about plotting figures, if you want it




%% Bias wrapper functions

% simulate data with a radial bias towards screen centre
params2 = [.1 .1 20 .2]; % 4th parameter is mu = radial bias
model = WithRadialBias2D(SwapModel2D); % add radial bias to SwapModel
simulatedData.errors = SampleFromModel2D(model, params2,...
        [1 numTrials], simulatedData); % simulate responses
fit3 = MLE(simulatedData, model); % fit


%% 2AFC fitting

model = TwoAFC2D(StandardMixtureModel2D); 
% use TwoAFC2D wrapper on model without misbinding (for speed)

dimensions = [1366; 768]; % x and y dimensions in units (pixels here)
simulatedData = GenerateDisplays2D(numTrials, ...
                    itemsPerTrial,2,dimensions); % mode = 2 for 2AFC
simulatedData.afcCorrect = SampleFromModel2D(model, [.1 20],...
                    [1 numTrials], simulatedData); 
					% simulate which 2AFC choices were correct – note only 2 params

fit2AFC = MLE(simulatedData, model); % ensure model has TwoAFC2D() wrapper


%% response sampling

% instead of assuming a uniform guessing distribution, you can build a
% response probability function across a person's responses and sample from
% this to get a response distribution, using the WithResponseSampling2D()
% wrapper function:

model = WithResponseSampling2D(SwapModel2D()); % no extra parameter

simulatedData = GenerateDisplays2D(numTrials, itemsPerTrial); % make stimuli
simulatedData.errors = SampleFromModel2D(model, params,...
                    [1 numTrials], simulatedData); % simulate responses

% you need the response coordinates on each trial
simulatedData.resps = simulatedData.errors - simulatedData.targets;

% you also need to store this in a cell so that it cannot be split (for
% example by SplitDataByCondition() )
simulatedData.respPdf = {simulatedData.resps};                
                
fit4 = MLE(simulatedData, model); % ensure model has TwoAFC2D() wrapper


%% stimulus constraints
% you can set constraints on where the stimuli may appear during learning,
% such as minimum distances from edges or other stimuli

% distances between [stimuli, stimuli & edges, stimuli & centre]
minDists = [88 88 29]; % pix equiv of [3 3 1] vis deg at 40cm
model = SwapModel2D(); 

% simulate with those constraints
simulatedData = GenerateDisplays2D(numTrials, itemsPerTrial,...
                        1, dimensions, minDists);
simulatedData.errors = SampleFromModel2D(model, params,...
                    [1 numTrials], simulatedData); % simulate responses

fit5 = MLE(simulatedData, model);

%% resample or bound
% the default is edge-resampling, but can be switched to edge-constraining

model = SwapModel2D(1); % set boundary to 1 when calling the model

simulatedData = GenerateDisplays2D(numTrials, itemsPerTrial,...
                        1, dimensions); % simulate as normal
simulatedData.errors = SampleFromModel2D(model, params,...
                    [1 numTrials], simulatedData); % simulate responses

fit6 = MLE(simulatedData, model); % fit as normal


