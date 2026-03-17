clear; clc; cla; clf

[sampleStruct, probStruct, Comments] = scfread('sample.scf');
% figure
% hold on
% plot(sampleStruct.A);
% plot(sampleStruct.C);
% plot(sampleStruct.G);
% plot(sampleStruct.T);

% legend('A','C','G','T');

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
