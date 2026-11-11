S1 = load("D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\Noisy\Matlab\Hela_ROI_07_30_7395_5114_161\ROI_07_025_Nuclei.mat");
noisy_nuclei = S1.Hela_nuclei;
noisy_bg = S1.Hela_background;
noisy_cell = load("D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\Noisy\Matlab\Hela_ROI_07_30_7395_5114_161\ROI_07_025_Cell.mat");

S2 = load("D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\Noisy\Matlab\Hela_ROI_07_30_7395_5114_161\ROI_07_000_Nuclei.mat");
healthy_nuclei = S2.Hela_nuclei;
healthy_bg = S2.Hela_background;
healthy_cell = load("D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\Noisy\Matlab\Hela_ROI_07_30_7395_5114_161\ROI_07_000_Cell.mat");

slice=100;

figure;
imshow(healthy_nuclei(:,:,slice));
title('Healthy');

figure;
imshow(noisy_nuclei(:,:,slice));
title('Noisy');