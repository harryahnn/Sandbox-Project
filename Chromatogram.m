clear; clc; cla; clf

[sampleStruct, probStruct, Comments] = scfread('sample.scf');
% figure
% hold on
% plot(sampleStruct.A);
% plot(sampleStruct.C);
% plot(sampleStruct.G);
% plot(sampleStruct.T);

% legend('A','C','G','T');

%%MAIN PORTION

newSampleStructA = removeNoise(sampleStruct.A);
newSampleStructC = removeNoise(sampleStruct.C);
newSampleStructG = removeNoise(sampleStruct.G);
newSampleStructT = removeNoise(sampleStruct.T);
figure;
hold on;
plot(newSampleStructA);
plot(newSampleStructC);
plot(newSampleStructG);
plot(newSampleStructT);
legend('A (Filtered)', 'C (Filtered)', 'G (Filtered)', 'T (Filtered)');

[prob,chain] = calculateProb(newSampleStructA, newSampleStructC, newSampleStructG, newSampleStructT);
chainStr = strjoin(chain, '');

gcPercent = calculateGCContent(chainStr);

%% FUNCTIONS



function newSampleStruct = removeNoise(sampleData)
    startPoint = 1; %creates the start point to cut the data from
    trimThreshold = mean(sampleData) * 0.1; %creates a threshold for "valuable" data

    for i = 1:length(sampleData) 
        if sampleData(i) > trimThreshold %if data is above the threshold to be considered "usefUl"
            startPoint = i; %we then update the starting point of the data
            break; %stop iterating because we have the starting point we want
        end
    end
    
    sampleData = sampleData(startPoint:end); %cut off unnecessary data at the beginning

    x = 1:length(sampleData); %creates our x values
    x = x'; %turns it into a row array

    thresholdValue = mean(sampleData) * 0.7; %creates a threshold for what values are considered "noise"
    noiseIndices = find(sampleData < thresholdValue); %finds the indices on each sampleData to find when the value is just noise
    noiseValues = sampleData(noiseIndices); %finds the actual value of the signal strength for that base pair

    p = polyfit(noiseIndices, noiseValues, 5); %creates parameters for the noise values of the sample (excludes all useful values)
    noiseThreshold = polyval(p, x); %creates a rough curve of what the noise is using parameters

    newSampleStruct = sampleData - noiseThreshold; %subtracts the "noise" from the actual data

    newSampleStruct(newSampleStruct < 0) = 0; %any negative values turn positive
end

function [baseProb, baseChain] = calculateProb(structA, structC, structG, structT) %finds probabilities(baseProb) of base pairs as well as the sequence(baseChain)

    minLen = min([length(structA), length(structC), length(structG), length(structT)]); %ensures that every single data set for the signal strength of each base pair is equal length
    structA = structA(1:minLen);
    structC = structC(1:minLen);
    structG = structG(1:minLen);
    structT = structT(1:minLen);

    globalMean = mean([structA; structC; structG; structT]); %finds the mean signal strength among each basepair
    threshold = globalMean * 0.5; %creates a threshold value that filters out all peaks below a certain height

    [ampA, locA] = findpeaks(structA, 'MinPeakHeight', threshold, 'MinPeakDistance', 5); %finds the amplitude and the location of every peak for each base pair
    [ampC, locC] = findpeaks(structC, 'MinPeakHeight', threshold, 'MinPeakDistance', 5); %ensures ample distance between pairs to avoid constant repeats
    [ampG, locG] = findpeaks(structG, 'MinPeakHeight', threshold, 'MinPeakDistance', 5); %also ensures that signals not strong enough are not counted
    [ampT, locT] = findpeaks(structT, 'MinPeakHeight', threshold, 'MinPeakDistance', 5); %this is our countermeasure against dye bleeding between base pairs

    allBases = table; %creates the table allBases with locations, amplitudes, and the sequence
    allBases.locations  = [locA;  locC;  locG;  locT];
    allBases.amplitudes = [ampA;  ampC;  ampG;  ampT];
    allBases.basepairs  = [repmat('A', length(locA), 1); repmat('C', length(locC), 1); repmat('G', length(locG), 1); repmat('T', length(locT), 1)]; %will return a sequence that looks like 'A A A C C C C C C G G G G G G G T T T T T T'

    allBases = sortrows(allBases, 'locations', 'ascend'); %sorts by location, so based off of where the location of each peak is, to get the proper sequence

    n = height(allBases); %used for the creation of our outputs
    baseChain = cell(1, n); %creates a cell array for the output for our sequence
    baseProb  = zeros(n, 4); %creates an array with 4 columns for the probabilities for every base pair at every location

    for i = 1:n %for every peak found
        pos = allBases.locations(i); %we get the position of the base 

        intensities = [structA(pos), structC(pos), structG(pos), structT(pos)]; %creates array of intensities using the amplitudes of the peaks at each position 
        intensities(intensities < 0) = 0; %if the intensity value at the peak is less than zero, set it to zero (doesn't happen, but just a safety measure)
        totalIntensity = sum(intensities); %gets the sum of the intensities at that point
        if totalIntensity > 0
            baseProb(i, :) = intensities / totalIntensity; %the probability for each base pair is their intensity over total intensity, and it plots this for each peak
        end

        baseChain{i} = allBases.basepairs(i); %sets the cell array baseChain to the sequence that was in allBases.
    end
end

function gcPercent = calculateGCContent(dnaSequence) % This function calculates the GC-content percentage in a DNA sequence. 
%Inputs: string containing the DNA sequence 
%Outputs: The percentage of G and C content in the DNA sequence

    totalLength = length(dnaSequence); %total length of DNA sequence
    
    dnaSequence = upper(dnaSequence); % Convert to uppercase to standardize

    if totalLength == 0
        error('DNA sequence has no valid bases.'); %safety measure if our DNA sequence is not usable
    end
    
    % Count the occurrences of G and C bases
    gcCount = count(dnaSequence, 'G') + count(dnaSequence, 'C'); %finds the total GC content value
    gcPercent = (gcCount/totalLength) * 100; %finds the total GC content value percentage

    windowSize = 50; %every window that we calculate the GC content from is 50 base pairs - less accurate, but less calculation than if we did a window size of 20
    
    xWindowData = []; %initialize arrays for our plotted data
    yWindowData = [];

    for i = 1:10:totalLength - windowSize %essentially, looks at a range of 50 base pairs from the length, gets the GC content, and then moves 10 base pairs forward
        window = dnaSequence(i:i + windowSize - 1); %sets the bounds for the window
        wgcCount = count(window, 'G') + count(window,'C'); %finds the window GC count
        wgcPercent = (wgcCount/windowSize) * 100; %finds the window GC percent
        xWindowData = [xWindowData, round((2 * i + windowSize)/2)]; %gets the center of every window
        yWindowData = [yWindowData, wgcPercent]; %gets every window GC percentage
    end

    figure;
    plot(xWindowData, yWindowData, 'b-', 'LineWidth', 1.5); hold on; %plots the actual curve for all the GC counts at every window
    yline(gcPercent, 'r--', sprintf('Overall: %.1f%%', gcPercent), 'LineWidth', 1.2); %creates a horizontal line for the average GC percent
    yline(50, 'k:', '50%', 'LineWidth', 1); %creates a horizontal line right at 50%
    xlabel('Position in sequence (bp)'); 
    ylabel('GC content (%)');
    title('GC content Curve');
    legend('Window GC%', 'Overall GC%', '50% Reference Line');
    ylim([0 100]); %0-100%

end  

function mutations = findMutations(seq1, seq2)

    len = min(length(seq1), length(seq2));
    mutationIndex = [];
    healthyBase = [];
    patientBase = [];

    for i = 1:len
        if seq1(i) ~= seq2(i)
            mutationIndex(end + 1) = i;
            healthyBase(end + 1) = seq1(i);
            patientBase(end + 1) = seq2(i);
        end
    end
 
    mutations.index = mutationIndex;
    mutations.healthy = healthyBase;
    mutations.patient = patientBase;
end 
