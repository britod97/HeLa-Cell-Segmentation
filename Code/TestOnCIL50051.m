cd("E:\HeLa\HeLa-Cell-Segmentation-master\Code")

%%
outDir = 'E:\HeLa\Data\CIL50051\Cropped_Tiffs_2';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for slice = 30:199
    filename = sprintf('E:\\HeLa\\Data\\CIL50051\\Downsampled_Tiffs\\slice_%04d.tif', slice);
    Hela0 = imread(filename);

    cropped = Hela0(1201:3200, 250:2249);

    outName = sprintf('slice_%04d_cropped.tif', slice);
    imwrite(cropped, fullfile(outDir, outName));
end

%%

% nuclei = segmentNucleiHelaEM_3D('E:\\HeLa\\Data\\CIL50051\\Cropped_Tiffs');

%%
slices = [170];
n = numel(slices);
%figure
for i = 1:n
    slice = slices(i);

    filename = sprintf('E:\\HeLa\\Data\\CIL50051\\Downsampled_Tiffs\\slice_%04d.tif', slice);
    Hela0 = imread(filename);

    ROI = Hela0(1201:3200,250:2249);
    Hela = double(ROI(:,:,1));
    Hela_nuclei = segmentNucleiHelaEM(Hela);

    subplot(2,n,i)
    imagesc(ROI); axis image off
    title(sprintf('ROI %d', slice))

    subplot(2,n,i+n)
    imagesc(Hela_nuclei); axis image off
end

colormap(gray)

%%
baseDir = 'E:\HeLa\Data\CIL50051\Cropped_Tiffs_2\';
dir0 = dir(fullfile(baseDir, '*.tif*'));
numSlices = numel(dir0);

% Step 1: Find the slice with the largest detected nucleus
nucleusArea = zeros(numSlices, 1);
for k = 1:numSlices
    img = imread(fullfile(baseDir, dir0(k).name));
    try
        mask = segmentNucleiHelaEM(double(img), [], 4);
        nucleusArea(k) = sum(mask(:));
    catch
        nucleusArea(k) = 0;
    end
    fprintf('Slice %d: nucleus area = %d\n', k, nucleusArea(k));
end

% Step 2: Use the largest nucleus slice as the seed
[~, bestSlice] = max(nucleusArea);
fprintf('Best central slice: %d\n', bestSlice);

% Step 3: Run the 3D segmentation seeded from the best slice
[Hela_nuclei, Hela_background] = segmentNucleiHelaEM_3D(baseDir, bestSlice, 4);