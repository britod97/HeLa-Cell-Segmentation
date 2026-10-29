fpath = ['D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094/NewMitochondriaPredictions/' ...
    '    Tiffs/ROI_01_30_876_1665_81_z0001-mito-prediction.tif']

info = imfinfo(fpath);
numImages = numel(info);

mitochondria = zeros(info(1).Height, info(1).Width, numImages, 'uint16'); % or uint8 depending on data

for k = 1:numImages
    mitochondria(:,:,k) = imread(fpath, k);
end

save("D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\NewMitochondriaPredictions\" + ...
    "    Matlab/ROI_01_30_876_1665_81_" + ...
    "Mitochondria.mat", "mitochondria", '-v7.3');