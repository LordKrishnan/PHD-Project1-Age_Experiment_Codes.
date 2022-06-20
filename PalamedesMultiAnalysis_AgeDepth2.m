% Combine 2 or more files for analysis for the stereoacuity task
% clear the workspace to start to make sure there's nothing bad hanging
% around from previous runs.
clear
commandwindow


%Analysis Constants
InitSlope = 1; % init slope parameter. Change if can't get fit to work

InitGuessRate = 0.5; % initial parameter. only used if this parameter is free
InitLapseRate = 0;  % initial parameter. only used if this parameter is free 


% Get the files
% This calls uigetfile to load a group of data files you are interested in
% analysing. filename and pathname are cell arrays which store the
% folders (pathname) and files (filename) to the data files
[filename, pathname] = uigetfile('*.txt', 'Pick an M-file','MultiSelect','on');


% step through the filename array, and load the data from each file as a
% separate array called 'data'
for i = 1:1:numel(filename)
    %load the first file
    [header{i} data{i}] = hdrload(fullfile(pathname,filename{i}));
end

% Now, compile all of the data into one large array.
% First, initialize the compiled data array, as an empty variable.
compiledData = [];
% step through each cell of 'data' and add it to the compiled data array
for i = 1:1:numel(data);
    compiledData = [compiledData;data{i}];
end

% This just copies what we had to a new variable called 'data'. This was
% done to be compatible with some older code. But, is not necessary. Can
% just do all of the analysis on compiledData
data = [];
data = compiledData;
clear compiledData;

distanceCol = 2;
depthCol = 7;
correctCol = 11;

% shape = 'cumulative gaussian';
% prefs = batch('shape',shape,'n_intervals',2,'runs',1000);
% outputPrefs = batch('write_pa', 'pa', 'write_th', 'th');

options = optimset('fminsearch');
options.MaxFunEvals = 2000;
options.MaxIter = 2000;
options.TolFun = 1e-09;             %require higher precision on LL
options.TolX = 1e-09;               %require higher precision on parameter
                                    %estimates
options.Display = 'off';            %suppress fminsearch output

% PF = @PAL_CumulativeNormal;
PF = @PAL_CumulativeNormal;
PFI = @PAL_inverseCumulativeNormal;
% PF=@PAL_Logistic
% PF = @PAL_Weibull;
numSimulations = 200; % number of simulations for bootstrapped SE's


% you pass paramaters [alpha beta lambda theta] to the fitting engine.
% These are, [threshold slope <guessing rate> <lapse rate>]
% The paramsFree vector is where you decide what you want to freely vary in
% the fit. 0 = fixed, 1 = freely vary. So, here we let the fit engine
% search for the best slope and threshold paramaters, which leaving the
% guessing and lapse rate to be fixed. The actual values of these are set
% later in the code.

%paramsValues = [30 5 0.5 0];
paramsFree = [1 1 0 1];
%paramsFree = [1 1 1 1];


options = optimset('fminsearch');   %type help optimset
options.TolFun = 1e-09;     %increase required precision on LL
options.Display = 'off';    %suppress fminsearch messages
lapseLimits = [0 1];        %limit range for lambda
                            %(will be ignored here since lambda is not a
                            %free parameter)
                            
                            % extract each distance from the data file to analyse
% the functions separately. See >>help unique for an explanation of how
% this works. Also, see my matlab tutorial that I sent you.
distances = unique(data(:,distanceCol));

for d = 1:1:numel(distances)

    % We need to search through the data array and find all of the trials
    % where the target distance was  == to the distance we're interested in
    %
    % Initialize an indexing array. Ask yourself, why is it being
    % initialized to empty on each iteration?
    dataIND = [];
    % find all of the entries in the data array equal to the distance we're
    % analysing this time through the loop. Store these indices in variable
    % dataIND
    dataIND = find(data(:,distanceCol)==distances(d));
    % Extract the data that is at each of those indices.
    distanceData{d} = data(dataIND,:);

    % Now, we take the data we've got for distance(d) and figure out what
    % levels of the stimulus the observer saw.
    %
    % This does the exact same thing as the stuff above, except on the
    % corrugation frequency data.
    depthLevels{d} = unique(distanceData{d}(:,depthCol));
    % step through each level of the stimulus
    for st = 1:1:numel(depthLevels{d});
        % initialize an empty variable for the INDEXING
        IND = [];
        % find the trials that we want, in this case its all the unique
        % depth levels 
        IND = find(distanceData{d}(:,depthCol)==depthLevels{d}(st));

        % how many trials did they get at that stimulus value?
        numtrials = numel(IND);
        
        % how many of those trials did they get correct?
        numcorrect = sum(distanceData{d}(IND,correctCol));
        % compile all of this data into an array called psychometric data.
        psychometricData{d}(st,:) = [depthLevels{d}(st) numcorrect numtrials];
    end
    
    % Now, we have all of the psychometric data compiled. Lets do some
    % function fitting.
    

    % Set the initial guesses of the psychometric function. This is just an
    % initial first guess, and the estimates of the threshold and slope are
    % widly wrong. They are seed values for the fminsearch. See >>help
    % fminsearch for an explanation of what happens next
    paramsValues{d} = [mean(psychometricData{d}(:,1)) InitSlope InitGuessRate InitLapseRate];
    paramsInput{d} = paramsValues{d};    
    % create a variable that ranges from the min and max of the stimulus
    % levels for plotting the psychometric function
    fineGrainedStimList = min(psychometricData{d}(:,1))-0.2:0.01:max(psychometricData{d}(:,1))+0.2;

    % fit the psychometric function. See >>help PAL_PFML_Fit for details of
    % the input parameters and the output values 
    [paramsValues{d} LL exitFlag output] = PAL_PFML_Fit(psychometricData{d}(:,1),...
        psychometricData{d}(:,2),...
        psychometricData{d}(:,3),...
        paramsValues{d},paramsFree,PF,'searchOptions',options,...
        'lapseLimits',lapseLimits);
    LapseLog(d) = paramsValues{d}(4);
    
    
    % open a figure window
%     figure
    
    % create a simple plot which has elements scaled by the number of
    % trials at some stimulus level.
    for i = 1:1:numel(psychometricData{d}(:,1))
        
        plot(psychometricData{d}(i,1),psychometricData{d}(i,2)./psychometricData{d}(i,3),'b.','markersize',psychometricData{d}(i,3).*1.1);
        hold on
        % compute values along the psychometric function given the best fit
        % parameters and the finegrained stimulus value vector. Plot the
        % psychometric function
        ProportionModel = PF(paramsValues{d},fineGrainedStimList);%,'g-','linewidth',4);
        subplot(3,1,d),plot(fineGrainedStimList,ProportionModel,'r-','linewidth',2);
        ylim([0.4 1.1]);
        xlim([min(psychometricData{d}(:,1))-.2 max(psychometricData{d}(:,1))+.2])
        xlabel('Stimulus Depth (seconds of arc)');
        ylabel('Proportion correct')
    end
    
    
    
    
    % We've extracted threshold and plotted the psychometric function, so
    % now, bootstrap some error estimates.
    % See >>help PAL_PFML_BootstrapNonParametric for details about what's
    % going on. Also, see >>help PAL_PFML_BootstrapParametric to see the
    %     % difference between a non-parametric and parametric bootstrap
    %     fprintf('Determining standard errors for distance %i of %i......\n')
    [SD paramsSim{d} LLSim converged] = PAL_PFML_BootstrapNonParametric(...
        psychometricData{d}(:,1),...
        psychometricData{d}(:,2),...
        psychometricData{d}(:,3),...
        paramsValues{d},paramsFree,numSimulations,PF);

    % compute the threshold by computing the inverse cumulative normal for
    % the 0.75 proportion correct section
    threshold(d) = PFI(paramsValues{d},0.90);
    %           threshold(d) = paramsValues(1)
    thresholdSD(d) = SD(1);

    % This computes the slope of the function at the 0.5 correct point
    slope(d) = paramsValues(d);
    slopeSD(d) = SD(d);


end


% Now we've got everyFixData we need for each distance run. So, compile it
% into an array called thresholdData and slopeData

% thresholdData now stores the threshold data in the first row and the
% estimated standard deviation in the second row.
thresholdData(1,:) = threshold;
thresholdData(2,:) = thresholdSD


% There is now an estimate of the slope of the function as well, called
% slopeData. Stores the slope of the psychometric function as well as it's
% estimated standard deviation from the bootstrap.
slopeData(1,:) = slope;
% slopeData(2,:) = slopeSD;
