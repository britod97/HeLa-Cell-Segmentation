% Input and output directories
inputDir = 'E:\HeLa\Data\EMPIAR-10094\Noisy\Predictions\Hela_ROI_08_30_6001_6460_221';
outputDir = 'E:\HeLa\Data\EMPIAR-10094\Noisy\Matlab\Hela_ROI_08_30_6001_6460_221';

% Create output directory if it doesn't exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Find all TIFF files
tifFiles = dir(fullfile(inputDir, '*.tif'));

for i = 1:length(tifFiles)

    % Full path to current TIFF
    fpath = fullfile(inputDir, tifFiles(i).name);

    % Get TIFF information
    info = imfinfo(fpath);
    numImages = numel(info);

    % Preallocate using the TIFF datatype
    switch info(1).BitDepth
        case 8
            dtype = 'uint8';
        case 16
            dtype = 'uint16';
        otherwise
            dtype = 'double';
    end

    volume = zeros(info(1).Height, info(1).Width, numImages, dtype);

    % Read each slice
    for k = 1:numImages
        volume(:,:,k) = imread(fpath, k);
    end

    % Construct output filename
    [~, baseName, ~] = fileparts(tifFiles(i).name);
    savePath = fullfile(outputDir, [baseName '.mat']);

    % Save variable
    save(savePath, 'volume', '-v7.3');

    fprintf('Converted %s -> %s\n', tifFiles(i).name, [baseName '.mat']);

end

disp('Done.');