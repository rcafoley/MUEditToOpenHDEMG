%Code for converting MUEdit exports to OTB style exports for use in
%openHDEMG and analyze HDEMG decompositions.
% 
%In openHDEMG open the new file as an OTB format and ensure
%you select the same extension factor that was chosen in MUEdit as this
%will not be hardcoded anywhere. 

%Note that this is designed for a very specific export. you need to change
%the descriptions dection to match your export or the torque data to pull
%your correct aux channel info. Further, this only imports the first
%channel of aux data as 'aquired force' if you wish to use some other
%channel as the reference signal then that must be changed to select the
%proper column of data from signal.auxillary on line 83

%Last edit April 4 2024

%Copyright Ryan C. A. Foley 2024
%Github: rcafoley
%Twitter: foleyneuromech
%LinkedIn: ryancafoley

%% Clear all variables from the workspace, clear the Command Window, and print the current directory
clear all; % Clears all variables from the workspace
clc;       % Clears the Command Window
disp(['Current working directory: ', pwd]); % Prints the current directory

%Variables to change to accomodate different experimental setups
%TypeofArray(s)
emgArray = 'GR04MM1305'
%muscle
muscle = 'Tibialis Anterior'
%MVC value & offset
mvc = 0.466 %in mV
offset = 0.13 %inmV

% List all .mat files in the current directory
files = dir('*.mat');

for i = 1:length(files)
    originalFileName = files(i).name;
    fprintf('Processing %s\n', originalFileName);
    
    % Load the .mat file
    loadedData = load(originalFileName);
    % Assume 'edition' and 'signal' are variables within 'loadedData'
    % Adjust the following lines if they're nested within another structure
    edition = loadedData.edition; % Or the appropriate field if nested differently
    signal = loadedData.signal;   % Adjust based on your data structure
    

    %% Data Conversion here

    %define the number of units from the decomp - this may need to the "clean"
    %versions if a unit is deleted
    
    numRawChannels = size(signal.data, 1);
    numMUs = size(edition.Dischargetimes, 2);

    %define the total channels to be created Array size (monopolar + 1
    %force + number of MUs x2 decomp and source)

    totalChannels = numRawChannels + 1 + numMUs + numMUs; 

    %FileType
    OTBFile = 'unknown';
    %SamplingRate
    SamplingFrequency = signal.fsamp;

    %Time Vector
    % Initialize the cell array
    Time = cell(1,1); % Creates a 1x1 cell array
    % Insert the array 'edition.time' into the {1,1} location of the cell array 'Time'
    Time{1,1} = edition.time';

    %Convert the HDEMG Array Data Arrays
    %first 64 are the raw data from the grid
    Data = cell(1,1); % Creates a 1x1 cell array
    Data{1,1} = signal.data';

    % Extract the specific torque data from 'signal.auxillary' - ours is in
    % row 1 but you may have more or want a different row for other aux or
    % the target path
   torqueColumnData = signal.auxiliary(10, :).'; % Transpose it to match the orientation (rows, not columns)
    
    % Check to ensure dimensions match for concatenation
    if size(Data{1,1}, 1) == length(torqueColumnData)
        % Add the new column to the array in 'Data{1,1}'
        Data{1,1} = [Data{1,1}, torqueColumnData];
    else
        disp('Error: Dimension mismatch. Cannot concatenate torque to the new column.');
    end

    %Extract the decomposition aka discharge times
    maxSamples = length(edition.time);

    % Initialize the array to hold all converted data. Rows represent samples, columns represent scenarios.
    allDecompositionConvertedData = zeros(maxSamples, numMUs);
    
    % Loop through each scenario and fill in the allConvertedData array
    for i = 1:numMUs
        % Extract the event data for the current scenario
        spikeData = edition.Dischargetimes{1, i};
        
        % Create a temporary array of zeros for the current scenario
        otbSpikeData = zeros(maxSamples, 1);
        
        % Mark the events in the tempData
        otbSpikeData(spikeData) = 1;
        
        % Assign the converted data to the corresponding column in allConvertedData
        allDecompositionConvertedData(:, i) = otbSpikeData;
    end

        % Check to ensure dimensions match for concatenation
    if size(Data{1,1}, 1) == length(allDecompositionConvertedData)
        % Add the new column to the array in 'Data{1,1}'
        Data{1,1} = [Data{1,1}, allDecompositionConvertedData];
    else
        disp('Error: Dimension mismatch. Cannot concatenate spike decomposition to the new column.');
    end

        %Extract the source data for the decomposition aka individual optimised
    %EMG channels for that MU **** NOTE: this may have to be changed to
    %PulseTrainCLEAN as that may be the edited data - same for decomp
    sourceEMGColumnData = edition.Pulsetrain{1,1}.'; % Transpose it to match the orientation (rows, not columns)
    
    % Check to ensure dimensions match for concatenation
    if size(Data{1,1}, 1) == length(sourceEMGColumnData)
        % Add the new column to the array in 'Data{1,1}'
        Data{1,1} = [Data{1,1}, sourceEMGColumnData];
    else
        disp('Error: Dimension mismatch. Cannot concatenate sourceEMG to the new column.');
    end

    %Reconstruct the needed Descriptions for each channel
    %these cannot be extracted from the MUedit export easily so we will
    %reconstruct them. If your setup is different you will need to adjust
    %these to have the coorect info. this can be found 

    % Initialize the Descriptions cell array
    Description = cell(totalChannels, 1); % Adjust the size based on the total number of descriptions
    
    % Fill the array with descriptions using dynamic components
    for i = 1:numRawChannels % For the repetitive part with numbers 1 to 64
        Description{i} = sprintf('%s - MULTIPLE IN 1 (Channel 1->64) - %s (%d)[mV]', muscle, emgArray, i);
    end
    
    % Add unique descriptions after the loop with dynamic MVC and Offset values
    Description{numRawChannels+1} = sprintf('acquired data MVC=%.3f Offset=%.2f[mV]', mvc, offset);
    
    % For Decomposition and Source for decomposition descriptions using dynamic components
    for i = 1:numMUs
        Description{(numRawChannels+1) + i} = sprintf('Decomposition of %s - MULTIPLE IN 1 (Channel 1->64) - %s (%d)[a.u]', muscle, emgArray, i);
        Description{(numRawChannels+1+numMUs) + i} = sprintf('Source for decomposition of %s - MULTIPLE IN 1 (Channel 1->64) - %s (%d)[a.u]', muscle, emgArray, i);
    end
    
    % Display the Descriptions cell array to verify
    %disp(Description);

    %% Output here
    %% Output here
    % Save the original data with a new name (_MUEditExport appended)
    originalFileNewName = [originalFileName(1:end-4), '_MUEditExport.mat'];
    save(originalFileNewName, '-struct', 'loadedData');
    
    % Save the variables into the new MAT file
    newFileName = [originalFileName(1:end-4), '_ConvertedToOTBExportFormat.mat'];
    save(newFileName, 'Time', 'Data', 'Description', 'SamplingFrequency', 'OTBFile');
end