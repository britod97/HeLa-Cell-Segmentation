% Load cell membrane
membrane = load("D:\GitHub Repos\HeLa_Cell_Data\CIL50051\Matlab\ROI_01_Cell.mat");

baseDir = 'D:\GitHub Repos\HeLa_Cell_Data\CIL50051\Bigger_Cropped_Tiffs';
dir0 = dir(fullfile(baseDir, '*.tif*'));
numSlices = numel(dir0);

% Read first slice to get dimensions
info = imfinfo(fullfile(baseDir, dir0(1).name));
Hela_nuclei = false(info.Height, info.Width, numSlices);

% se = strel("disk",5);

for k = 1:numSlices
    img = double(imread(fullfile(baseDir, dir0(k).name)));
    % membraneSlice = imdilate(membrane.mask(:,:,k), se);
    membraneSlice = membrane.mask(:,:,k);
    % imgMasked = img.*membraneSlice;
    % segmentation = segmentNucleiHelaEM(imgMasked, [], 4);
    segmentation = segmentNucleiHelaEM(img, [], 5);
    Hela_nuclei(:,:,k) = logical(segmentation);
    fprintf('Processed slice %d/%d\n', k, numSlices);
end

% Apply the same 3D median filter the original code uses
Hela_nuclei = medfilt3(Hela_nuclei, [3 3 13]);

volshow(Hela_nuclei);
% Save result
% save('ROI_01_Nuclei_masked.mat', 'Hela_nuclei', '-v7.3');
%%
% Hela_nuclei = load("nuclei_segmentation_cefas.mat")
% arr = Hela_nuclei.Hela_nuclei;
% 
% volshow(arr)
%%
% factor = 1;  % 2000 -> 500 in XY
% arr_small = arr(1:factor:end, 1:factor:end, :);
% 
% % Smooth slightly to get a nicer surface (optional but helps)
% % arr_small = smooth3(double(arr_small), 'gaussian', 3);
% 
% volshow(arr)
% % figure
% % p = patch(isosurface(arr_small, 0.5));
% % %isonormals(arr, p)
% % set(p, 'FaceColor', 'red', 'EdgeColor', 'none')
% 
% daspect([1 1 1])
% view(3)
% axis tight
% camlight
% lighting gouraud
% title('3D visualization of nuclei mask')
