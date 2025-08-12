%Code for converting MUEdit exports to OTB style exports for use in
%openHDEMG and analyze HDEMG decompositions.
% Now supports both 1 and 2 electrode array configurations
% 
%In openHDEMG open the new file as an OTB format and ensure
%you select the same extension factor that was chosen in MUEdit as this
%will not be hardcoded anywhere. 

%Note that this is designed for a very specific export. you need to change
%the descriptions section to match your export or the torque data to pull
%your correct aux channel info. 

%Last edit: Modified for 2-electrode support based on actual data structure
%Copyright Ryan C. A. Foley 2024
%Modified for multi-electrode support
%Github: rcafoley
%Twitter: foleyneuromech
%LinkedIn: ryancafoley

%% Clear all variables from the workspace, clear the Command Window, and print the current directory
clear all; % Clears all variables from the workspace
clc;       % Clears the Command Window
disp(['Current working directory: ', pwd]); % Prints the current directory

%% Configuration for different experimental setups
% MODIFY THESE FOR YOUR SETUP
% For single electrode setup, use same values for both or leave second empty
emgArrays = {'GR04MM1305', 'GR04MM1306'}; % Add second array ID here for 2-electrode setup
muscles = {'FCR', 'ECR'}; % Add second muscle name if different

% MVC values & offsets (add second values if different for each array)
mvc_values = [0.466, 0.500]; % in mV
offset_values = [0.13, 0.15]; % in mV

% Auxiliary channel to use for torque/force (1-indexed)
auxChannel = 10; % Change this to match your auxiliary channel (was row 10 in original)

% Channels per array (typically 64)
channelsPerArray = 64;

%% Process files
% List all .mat files in the current directory
files = dir('*.mat');

for fileIdx = 1:length(files)
    originalFileName = files(fileIdx).name;
    fprintf('\n========================================\n');
    fprintf('Processing %s\n', originalFileName);
    fprintf('========================================\n');
    
    % Load the .mat file
    loadedData = load(originalFileName);
    edition = loadedData.edition;
    signal = loadedData.signal;
    
    %% Detect number of electrode arrays
    numElectrodes = size(edition.Dischargetimes, 1);
    fprintf('Detected %d electrode array(s)\n', numElectrodes);
    
    % Get time and sampling info
    maxSamples = length(edition.time);
    SamplingFrequency = signal.fsamp;
    
    %% Initialize storage
    allDescriptions = {};
    
    %% Process Raw EMG Data
    fprintf('Processing raw EMG data...\n');
    rawData = signal.data'; % Transpose to get samples x channels
    numTotalRawChannels = size(rawData, 2);
    
    % Create descriptions for raw channels
    if numElectrodes == 1
        % Single electrode
        for i = 1:numTotalRawChannels
            description = sprintf('%s - MULTIPLE IN 1 (Channel 1->%d) - %s (%d)[mV]', ...
                muscles{1}, channelsPerArray, emgArrays{1}, i);
            allDescriptions{end + 1} = description;
        end
    else
        % Two electrodes - first 64 channels from array 1, next 64 from array 2
        for elecIdx = 1:numElectrodes
            for i = 1:channelsPerArray
                globalChannelNum = (elecIdx - 1) * channelsPerArray + i;
                description = sprintf('%s - MULTIPLE IN %d (Channel 1->%d) - %s (%d)[mV]', ...
                    muscles{elecIdx}, elecIdx, channelsPerArray, emgArrays{elecIdx}, i);
                allDescriptions{end + 1} = description;
            end
        end
    end
    
    %% Process Torque/Force Data
    fprintf('Adding torque/force data from auxiliary channel %d...\n', auxChannel);
    torqueData = signal.auxiliary(auxChannel, :)';
    
    % Add torque description
    if numElectrodes == 1
        torqueDesc = sprintf('acquired data MVC=%.3f Offset=%.2f[mV]', ...
            mvc_values(1), offset_values(1));
    else
        torqueDesc = sprintf('acquired data MVC=[%.3f,%.3f] Offset=[%.2f,%.2f][mV]', ...
            mvc_values(1), mvc_values(2), offset_values(1), offset_values(2));
    end
    allDescriptions{end + 1} = torqueDesc;
    
    %% Process Decomposition Data
    fprintf('Processing decomposition data...\n');
    allDecompData = [];
    actualMUsPerElectrode = zeros(1, numElectrodes);
    
    for elecIdx = 1:numElectrodes
        % Count actual MUs (non-empty cells) for this electrode
        muCount = 0;
        for j = 1:size(edition.Dischargetimes, 2)
            if ~isempty(edition.Dischargetimes{elecIdx, j})
                muCount = muCount + 1;
            end
        end
        actualMUsPerElectrode(elecIdx) = muCount;
        fprintf('  Electrode %d: %d MUs detected\n', elecIdx, muCount);
        
        % Create decomposition data for this electrode
        decompData = zeros(maxSamples, muCount);
        muIdx = 0;
        
        for j = 1:size(edition.Dischargetimes, 2)
            if ~isempty(edition.Dischargetimes{elecIdx, j})
                muIdx = muIdx + 1;
                spikeData = edition.Dischargetimes{elecIdx, j};
                otbSpikeData = zeros(maxSamples, 1);
                otbSpikeData(spikeData) = 1;
                decompData(:, muIdx) = otbSpikeData;
            end
        end
        
        allDecompData = [allDecompData, decompData];
    end
    
    % Add decomposition descriptions
    muGlobalCounter = 0;
    for elecIdx = 1:numElectrodes
        for muIdx = 1:actualMUsPerElectrode(elecIdx)
            muGlobalCounter = muGlobalCounter + 1;
            description = sprintf('Decomposition of %s - MULTIPLE IN %d (Channel 1->%d) - %s (MU%d)[a.u]', ...
                muscles{min(elecIdx, length(muscles))}, elecIdx, channelsPerArray, ...
                emgArrays{min(elecIdx, length(emgArrays))}, muGlobalCounter);
            allDescriptions{end + 1} = description;
        end
    end
    
    %% Process Source EMG Data
    fprintf('Processing source EMG data...\n');
    allSourceData = [];
    
    if numElectrodes == 1
        % Single electrode - Pulsetrain{1,1} contains all MU sources
        sourceMatrix = edition.Pulsetrain{1, 1}';
        % Only take the MUs that have decomposition data
        allSourceData = sourceMatrix(:, 1:actualMUsPerElectrode(1));
    else
        % Two electrodes - Pulsetrain{1,1} and Pulsetrain{1,2} contain sources
        for elecIdx = 1:numElectrodes
            sourceMatrix = edition.Pulsetrain{1, elecIdx}';
            % Only take the MUs that have decomposition data
            sourceData = sourceMatrix(:, 1:actualMUsPerElectrode(elecIdx));
            allSourceData = [allSourceData, sourceData];
        end
    end
    
    % Add source descriptions
    muGlobalCounter = 0;
    for elecIdx = 1:numElectrodes
        for muIdx = 1:actualMUsPerElectrode(elecIdx)
            muGlobalCounter = muGlobalCounter + 1;
            description = sprintf('Source for decomposition of %s - MULTIPLE IN %d (Channel 1->%d) - %s (MU%d)[a.u]', ...
                muscles{min(elecIdx, length(muscles))}, elecIdx, channelsPerArray, ...
                emgArrays{min(elecIdx, length(emgArrays))}, muGlobalCounter);
            allDescriptions{end + 1} = description;
        end
    end
    
    %% Combine all data
    Data = cell(1, 1);
    Data{1, 1} = [rawData, torqueData, allDecompData, allSourceData];
    
    % Convert descriptions to column cell array
    Description = allDescriptions';
    
    %% Prepare other variables
    OTBFile = 'unknown';
    
    % Time Vector
    Time = cell(1, 1);
    Time{1, 1} = edition.time';
    
    %% Summary
    totalMUs = sum(actualMUsPerElectrode);
    fprintf('\n--- Conversion Summary ---\n');
    fprintf('Total channels in output: %d\n', size(Data{1,1}, 2));
    fprintf('  - Raw EMG channels: %d\n', numTotalRawChannels);
    fprintf('  - Force/Torque channels: 1\n');
    fprintf('  - Decomposition channels: %d\n', totalMUs);
    fprintf('  - Source EMG channels: %d\n', totalMUs);
    if numElectrodes == 2
        fprintf('  - MUs from electrode 1: %d\n', actualMUsPerElectrode(1));
        fprintf('  - MUs from electrode 2: %d\n', actualMUsPerElectrode(2));
    end
    
    %% Verify dimensions
    if size(Data{1,1}, 2) ~= length(Description)
        warning('Dimension mismatch! Data has %d columns but %d descriptions', ...
            size(Data{1,1}, 2), length(Description));
    end
    
    %% Output
    % Save the original data with a new name
    originalFileNewName = [originalFileName(1:end-4), '_MUEditExport.mat'];
    save(originalFileNewName, '-struct', 'loadedData');
    fprintf('\nOriginal data saved as: %s\n', originalFileNewName);
    
    % Save the converted data
    newFileName = [originalFileName(1:end-4), '_ConvertedToOTBExportFormat.mat'];
    save(newFileName, 'Time', 'Data', 'Description', 'SamplingFrequency', 'OTBFile');
    fprintf('Converted data saved as: %s\n', newFileName);
end

fprintf('\n========================================\n');
fprintf('All files processed successfully!\n');
fprintf('========================================\n');