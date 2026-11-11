rootDir      = 'D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\Noisy\Tiffs';
segmentsDir  = 'D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\Noisy\Matlab';

% Keep only real folders (drop '.' and '..' and any stray files)
rois = dir(rootDir);
rois = rois([rois.isdir] & ~ismember({rois.name}, {'.','..'}));

% for i = 1:numel(rois)
for i = 1:1
    ROI_name = rois(i).name;                         % e.g. Hela_ROI_07_30_7395_5114_161
    ROI_path = fullfile(rois(i).folder, ROI_name);
    ROI_id   = regexp(ROI_name, 'ROI_\d+', 'match', 'once');   % e.g. ROI_07
    display(ROI_path)
    noiseLevels = dir(ROI_path);
    noiseLevels = noiseLevels([noiseLevels.isdir] & ~ismember({noiseLevels.name}, {'.','..'}));

    outdir = fullfile(segmentsDir, ROI_name);
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % for j = 1:numel(noiseLevels)
    for j = 1:2
        noise_name = noiseLevels(j).name;            % e.g. 000
        noise_path = fullfile(noiseLevels(j).folder, noise_name);
        display(noise_path)
        noise_id   = regexp(noise_name, '\d+', 'match', 'once');

        cellpath   = fullfile(outdir, sprintf('%s_%s_Cell.mat',   ROI_id, noise_id));
        nucleipath = fullfile(outdir, sprintf('%s_%s_Nuclei.mat', ROI_id, noise_id));

        if exist(cellpath, 'file') && exist(nucleipath, 'file')
            fprintf('Skipping %s / %s (already done)\n', ROI_name, noise_name);
            continue
        end

        fprintf('Processing %s / %s\n', ROI_name, noise_name);

        [Hela_nuclei, Hela_background] = segmentNucleiHelaEM_3D(noise_path);
        Hela_nuclei     = logical(Hela_nuclei);
        Hela_background = logical(Hela_background);

        % Save nucleus + background first, so they are kept even if the cell step fails
        save(nucleipath, 'Hela_nuclei', 'Hela_background', '-v7.3');

        Hela_cell = segmentCellHelaEM_3D(Hela_nuclei, Hela_background);
        Hela_cell = logical(Hela_cell);

        save(cellpath, 'Hela_cell', '-v7.3');
    end
end