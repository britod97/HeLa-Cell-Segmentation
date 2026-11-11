function segment_roi_cells(ROI_name, rootDir, segmentsDir)
    ROI_id   = regexp(ROI_name, 'ROI_\d+', 'match', 'once');
    ROI_path = fullfile(rootDir, ROI_name);
    outdir   = fullfile(segmentsDir, ROI_name);
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    levels = dir(ROI_path);
    levels = levels([levels.isdir] & ~ismember({levels.name}, {'.','..'}));

    for j = 1:numel(levels)
        noise_name = levels(j).name;
        noise_path = fullfile(levels(j).folder, noise_name);
        noise_id   = regexp(noise_name, '\d+', 'match', 'once');
        cellpath   = fullfile(outdir, sprintf('%s_%s_Cell.mat',   ROI_id, noise_id));
        nucleipath = fullfile(outdir, sprintf('%s_%s_Nuclei.mat', ROI_id, noise_id));

        if exist(cellpath, 'file') && exist(nucleipath, 'file')
            fprintf('Skipping %s / %s (already done)\n', ROI_name, noise_name);
            continue
        end

        try
            fprintf('Processing %s / %s\n', ROI_name, noise_name);
            [Hela_nuclei, Hela_background] = segmentNucleiHelaEM_3D(noise_path);
            Hela_nuclei     = logical(Hela_nuclei);
            Hela_background = logical(Hela_background);
            save(nucleipath, 'Hela_nuclei', 'Hela_background', '-v7.3');

            Hela_cell = logical(segmentCellHelaEM_3D(Hela_nuclei, Hela_background));
            save(cellpath, 'Hela_cell', '-v7.3');
        catch ME
            fprintf(2, 'FAILED %s / %s: %s\n', ROI_name, noise_name, ME.message);
        end
    end
end