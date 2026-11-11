function segment_noise_level(ROI_name, noise_name, rootDir, segmentsDir)
    ROI_id     = regexp(ROI_name, 'ROI_\d+', 'match', 'once');
    noise_path = fullfile(rootDir, ROI_name, noise_name);   % e.g. ...\Hela_ROI_08...\025
    noise_id   = regexp(noise_name, '\d+', 'match', 'once');
    outdir     = fullfile(segmentsDir, ROI_name);
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    cellpath   = fullfile(outdir, sprintf('%s_%s_Cell.mat',   ROI_id, noise_id));
    nucleipath = fullfile(outdir, sprintf('%s_%s_Nuclei.mat', ROI_id, noise_id));

    if exist(cellpath, 'file') && exist(nucleipath, 'file')
        fprintf('Skipping %s / %s (already done)\n', ROI_name, noise_name);
        return
    end

    fprintf('Processing %s / %s\n', ROI_name, noise_name);
    [Hela_nuclei, Hela_background] = segmentNucleiHelaEM_3D(noise_path);
    Hela_nuclei     = logical(Hela_nuclei);
    Hela_background = logical(Hela_background);
    save(nucleipath, 'Hela_nuclei', 'Hela_background', '-v7.3');

    Hela_cell = logical(segmentCellHelaEM_3D(Hela_nuclei, Hela_background));
    save(cellpath, 'Hela_cell', '-v7.3');
end