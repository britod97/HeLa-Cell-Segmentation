fpath = ['D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\NewMitochondriaPredictions\ROI_08_NOISY_1.5-mito-prediction.tif']

info = imfinfo(fpath);
numImages = numel(info);

mitochondria = zeros(info(1).Height, info(1).Width, numImages, 'uint16'); % or uint8 depending on data

for k = 1:numImages
    mitochondria(:,:,k) = imread(fpath, k);
end

save("D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\NewMitochondriaPredictions\Matlab\Noisy\ROI_08_NOISY_1.5" + ...
    "Mitochondria.mat", "mitochondria", '-v7.3');