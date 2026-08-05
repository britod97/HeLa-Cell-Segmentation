function cropTiffs(tiffdir, x, y, z_lo, z_hi, w, h, outdir)

    % Get all tif/tiff files
    files = dir(fullfile(tiffdir, '*.tif'));
    files = [files; dir(fullfile(tiffdir, '*.tiff'))];

    if isempty(files)
        error('No TIFF files found in directory: %s', tiffdir);
    end

    % Create output directory
    % outdir = fullfile(tiffdir, 'cropped');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % Crop each image
    for i = z_lo:z_hi
        
        infile = fullfile(tiffdir, files(i).name);

        % Read image
        I = imread(infile);

        % Crop
        Ic = histeq(imcrop(I, [x y w-1 h-1]));

        % Save
        outfile = fullfile(outdir, files(i).name);
        imwrite(Ic, outfile);

        fprintf('Cropped %s (%d/%d)\n', files(i).name, i, length(files));

    end

    fprintf('Finished cropping %d TIFFs.\n', length(files));

end