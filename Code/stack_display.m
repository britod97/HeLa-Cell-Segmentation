baseDir = 'D:\GitHub Repos\HeLa_Cell_Data\Data\Cropped_Tiffs_2';
dir0 = dir(fullfile(baseDir, '*.tif*'));
numSlices = numel(dir0);

% Read first slice to get dimensions
info = imfinfo(fullfile(baseDir, dir0(1).name));
Hela_nuclei = false(info.Height, info.Width, numSlices);

for k = 1:numSlices
    img = double(imread(fullfile(baseDir, dir0(k).name)));
    mask = segmentNucleiHelaEM(img, [], 4);
    Hela_nuclei(:,:,k) = logical(mask);
    fprintf('Processed slice %d/%d\n', k, numSlices);
end

% Apply the same 3D median filter the original code uses
Hela_nuclei = medfilt3(Hela_nuclei, [3 3 13]);

% Save result
save('nuclei_segmentation.mat', 'Hela_nuclei', '-v7.3');
%%
Hela_nuclei = load("nuclei_segmentation.mat")
arr = Hela_nuclei.Hela_nuclei;
%%
factor = 4;  % 2000 -> 500 in XY
arr_small = arr(1:factor:end, 1:factor:end, :);

% Smooth slightly to get a nicer surface (optional but helps)
arr_small = smooth3(double(arr_small), 'gaussian', 3);

figure
p = patch(isosurface(arr_small, 0.5));
%isonormals(arr, p)
set(p, 'FaceColor', 'red', 'EdgeColor', 'none')

daspect([1 1 1])
view(3)
axis tight
camlight
lighting gouraud
title('3D visualization of nuclei mask')
