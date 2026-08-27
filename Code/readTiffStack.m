fpath = 'D:/GitHub Repos/HeLa_Cell_Data/CIL50051/Mito/ROI_01_Mitochondria.tif'

info = imfinfo(fpath);
numImages = numel(info);

mitochondria = zeros(info(1).Height, info(1).Width, numImages, 'uint16'); % or uint8 depending on data

for k = 1:numImages
    mitochondria(:,:,k) = imread(fpath, k);
end

save("D:/GitHub Repos/HeLa_Cell_Data/CIL50051/Matlab/ROI_01_Mitochondria.mat", "mitochondria", '-v7.3');