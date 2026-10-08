function cropTiffsOctave(tiffdir, x, y, z_lo, z_hi, w, h, outdir)

    % Get all tif/tiff files, sorted by name so slice indices are stable
    files = dir(fullfile(tiffdir, '*.tif'));
    files = [files; dir(fullfile(tiffdir, '*.tiff'))];

    if isempty(files)
        error('No TIFF files found in directory: %s', tiffdir);
    end

    [~, order] = sort({files.name});
    files = files(order);

    if z_lo < 1 || z_hi > numel(files) || z_lo > z_hi
        error('z range [%d, %d] is invalid for %d files.', z_lo, z_hi, numel(files));
    end

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    for i = z_lo:z_hi

        infile = fullfile(tiffdir, files(i).name);
        I = imread(infile);

        % Crop by indexing (equivalent to imcrop(I, [x y w-1 h-1]))
        r2 = min(y + h - 1, size(I, 1));
        c2 = min(x + w - 1, size(I, 2));
        Ic = I(y:r2, x:c2, :);

        outfile = fullfile(outdir, files(i).name);
        imwrite(Ic, outfile);

        fprintf('Cropped %s (%d/%d)\n', files(i).name, i, numel(files));

    end

    fprintf('Finished cropping %d TIFFs.\n', z_hi - z_lo + 1);

end
