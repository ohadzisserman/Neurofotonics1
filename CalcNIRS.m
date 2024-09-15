function [dHbR, dHbO, fig] = CalcNIRS(dataFile, SDS, tissueType, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile, plotChannelIdx)

% Validate input parameters
if isfile(dataFile) && isfile(extinctionCoefficientsFile) && isfile(DPFperTissueFile) && isfile(relDPFfile) && isa(SDS, 'double') && isa(tissueType, 'char')
    % Set default value for plotChannelIdx if not provided
    if nargin < 7
        plotChannelIdx = [];
    end

    % Import intensity data
    rawData = load(dataFile);
    wavelengths = rawData.SD.Lambda;
    timeAxis = rawData.t;
    intensityMeasurements = rawData.d;

    % Extract extinction coefficients
    extinctionTable = readtable(extinctionCoefficientsFile);
    extCoeff1 = extinctionTable{extinctionTable.wavelength == wavelengths(1), {'HbO2', 'HHb'}};
    extCoeff2 = extinctionTable{extinctionTable.wavelength == wavelengths(2), {'HbO2', 'HHb'}};

    % Retrieve tissue-specific DPF
    tissueTable = readtable(DPFperTissueFile, 'Delimiter', '\t');
    tissueDPF = tissueTable{strcmp(tissueTable.Tissue, tissueType), 'DPF'};

    % Get wavelength-specific DPF coefficients
    dpfCoeffTable = readtable(relDPFfile);
    dpfCoeff1 = dpfCoeffTable{dpfCoeffTable.wavelength == wavelengths(1), 'relDPFcoeff'};
    dpfCoeff2 = dpfCoeffTable{dpfCoeffTable.wavelength == wavelengths(2), 'relDPFcoeff'};

    % Calculate effective DPF for each wavelength
    effectiveDPF1 = tissueDPF * dpfCoeff1;
    effectiveDPF2 = tissueDPF * dpfCoeff2;

    % Separate intensity data by wavelength
    intensity1 = intensityMeasurements(:, 1:20);
    intensity2 = intensityMeasurements(:, 21:40);

    % Compute change in optical density
    logIntensity1 = log(intensity1);
    logIntensity2 = log(intensity2);
    deltaOD1 = [-logIntensity1 + logIntensity1(1,:)];
    deltaOD2 = [-logIntensity2 + logIntensity2(1,:)];

    % Determine effective pathlength
    effectivePathLength1 = effectiveDPF1 * SDS;
    effectivePathLength2 = effectiveDPF2 * SDS;

    % Construct extinction coefficient matrix
    lambda1_HbO = extCoeff1(1);
    lambda1_HbR = extCoeff1(2);
    lambda2_HbO = extCoeff2(1);
    lambda2_HbR = extCoeff2(2);
    extinctionMatrix = [lambda1_HbO/effectivePathLength1, lambda1_HbR/effectivePathLength1; ...
                        lambda2_HbO/effectivePathLength2, lambda2_HbR/effectivePathLength2];

    % Invert extinction coefficient matrix
    invertedExtinctionMatrix = inv(extinctionMatrix);

    % Initialize concentration change matrices
    dHbO = zeros(size(deltaOD1));
    dHbR = zeros(size(deltaOD1));

    % Calculate concentration changes for each channel
    for channelIdx = 1:20
        deltaOD = [deltaOD1(:,channelIdx), deltaOD2(:,channelIdx)];
        concentrationChanges = invertedExtinctionMatrix * deltaOD';
        dHbO(:, channelIdx) = concentrationChanges(1,:);
        dHbR(:, channelIdx) = concentrationChanges(2,:);
    end

    % Generate plots if requested
    if ~isempty(plotChannelIdx)
        fig = figure;
        hold on;
        for channelIdx = plotChannelIdx
            subplot(length(plotChannelIdx), 1, find(plotChannelIdx == channelIdx));
            plot(timeAxis, dHbO(:,channelIdx), 'r', 'DisplayName', 'HbO');
            hold on;
            plot(timeAxis, dHbR(:,channelIdx), 'b', 'DisplayName', 'HbR');
            title(['Channel ' num2str(channelIdx)]);
            xlabel('Time (s)');
            ylabel('Concentration Change');
            legend show;
        end
        hold off;
    else
        fig = [];
    end

else
    error("Invalid input parameters. Please check file paths and data types.")
end
% Perform Fourier analysis on the first channel
firstChannelData = intensityMeasurements(:, 1);
dataLength = length(firstChannelData);
samplingRate = 1 / (timeAxis(2) - timeAxis(1));

% Compute FFT and shift it
fftResult = fft(firstChannelData);
fftResult = fftshift(fftResult);

% Create a properly centered frequency axis
if mod(dataLength, 2) == 0
    frequencyAxis = (-dataLength/2:dataLength/2-1) * (samplingRate/dataLength);
else
    frequencyAxis = (-(dataLength-1)/2:(dataLength-1)/2) * (samplingRate/dataLength);
end

% Visualize Fourier transform
figure;
plot(frequencyAxis, abs(fftResult));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform of First Channel');
xlim([-samplingRate/2, samplingRate/2]);  % Limit x-axis to Nyquist frequency range;

% Compute Signal-to-Noise Ratio (SNR)
assumedHeartRate = 1; % Hz, למשל. יש להתאים לפי הנתונים האמיתיים
noiseFrequencies = abs(frequencyAxis) > 2.5;

heartRateIndex = find(abs(frequencyAxis - assumedHeartRate) == min(abs(frequencyAxis - assumedHeartRate)), 1);
averageNoiseLevel = mean(abs(fftResult(noiseFrequencies)))
signalStrength = abs(fftResult(heartRateIndex))

SNR = signalStrength / averageNoiseLevel;
fprintf('Calculated Signal-to-Noise Ratio (SNR): %.2f\n', SNR);

end