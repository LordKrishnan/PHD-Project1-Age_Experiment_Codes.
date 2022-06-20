% Combine 2 or more files for analysis
clear
commandwindow
% Get the files
[filename, pathname] = uigetfile('*.txt', 'Pick an M-file','MultiSelect','on');



%load the first file
[header data] = hdrload(fullfile(pathname,filename));


% now compile all of the data into one big array
compiledData = [];
% for i = 1:1:numel(data);
compiledData = [compiledData;data];
% end

data = [];
data = compiledData;
clear compiledData;

% Instructions for how to read the imported text data. The distance of the
% target is in the third column. The corrugation frequency is in the 6th
% column, and whether the observer got the answer correct is in the 8th
% column
focalDistCol = 2;
distanceCol = 3;
stimSecCol = 7;
conflictCol = 8;
correctCol = 11;
Stair = 12;



% extract each distance from the data file to analyse
% the functions separately. See >>help unique for an explanation of how
% this works. Also, see my matlab tutorial that I sent you.
% conflicts = unique(data(:,conflictCol));
StairCases = unique(data(:,Stair));


% Everything is set up. Lets get fitting. Separate psychometric function
% per distance.
for c = 1:1:numel(StairCases)

    % We need to search through the data array and find all of the trials
    % where the target distance was  == to the distance we're interested in
    %
    % Initialize an indexing array. Ask yourself, why is it being
    % initialized to empty on each iteration?
    dataIND = [];
    % find all of the entries in the data array equal to the distance we're
    % analysing this time through the loop. Store these indices in variable
    % dataIND
    dataIND = find(data(:,Stair)==StairCases(c));
    % Extract the data that is at each of those indices.
    stairData{c} = data(dataIND,:);

    % Now, we take the data we've got for distance(c) and figure out what
    % levels of the stimulus the observer saw.
    %
    % This does the exact same thing as the stuff above, except on the
    % corrugation frequency data.
    arcSecLevels{c} = unique(stairData{c}(:,stimSecCol));
    % step through each level of the stimulus
    for st = 1:1:numel(arcSecLevels{c});
        % initialize an empty variable for the INDEXING
        IND = [];
        % find the trials that we want
        IND = find(stairData{c}(:,stimSecCol)==arcSecLevels{c}(st));
        % how many trials did they get at that stimulus value?
        numtrials = numel(IND);
        % how many of those trials did they get correct?
        numcorrect = sum(stairData{c}(IND,correctCol));
        % compile all of this data into an array called psychometric data.
        psychometricData{c}(st,:) = [arcSecLevels{c}(st) numcorrect numtrials];
    end
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the staircases for the data

for s = 1:1:c
    numtrials = [];
    numtrials{s} = size(stairData{s}(:,:));
    numtrials{s} = 1:1:numtrials{s}(:,1);
    
    figure(c+s)
    plot (numtrials{s}', stairData{s}(:,stimSecCol), 'b --o');
    
    hold off
end
    
