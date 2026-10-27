function cropAllROIs(tiffdir, outroot)
    % cropAllROIs  Crop a 2000x2000x300 region per ROI, using the ROI
    % coordinates encoded in the .mat filenames.
    %
    %   tiffdir : folder with the TIFF slices of the full volume
    %   outroot : parent folder; one subfolder per ROI is created inside

    matdir = 'E:\HeLa\Data\EMPIAR-10094\Matlab';
    W = 2000;  H = 2000;  D = 300;

    % ---- Find unique ROIs from the .mat filenames ----
    matfiles = dir(fullfile(matdir, 'Hela_ROI_*.mat'));
    if isempty(matfiles)
        error('No ROI .mat files found in: %s', matdir);
    end

    % Strip the trailing _Cell / _Invaginations / _Nuclei to get the ROI name
    roiNames = unique(regexprep({matfiles.name}, '_(Cell|Invaginations|Nuclei)\.mat$', ''));

    % ---- Volume dimensions (from the TIFF stack itself) ----
    files = [dir(fullfile(tiffdir, '*.tif')); dir(fullfile(tiffdir, '*.tiff'))];
    if isempty(files)
        error('No TIFF files found in directory: %s', tiffdir);
    end
    [~, idx] = sort({files.name});
    files = files(idx);
    nSlices = numel(files);                       % expected: 517

    info = imfinfo(fullfile(tiffdir, files(1).name));
    imgW = info(1).Width;
    imgH = info(1).Height;

    % ---- Process each ROI ----
    for k = 1:numel(roiNames)
        roi = roiNames{k};

        % Parse: ...<x>_<y>_<z> -> last three numeric tokens in the name
        tok = regexp(roi, '_(\d+)_(\d+)_(\d+)$', 'tokens', 'once');
        if isempty(tok)
            warning('Could not parse coordinates from "%s", skipping.', roi);
            continue
        end
        xc = str2double(tok{1});   % middle of x
        yc = str2double(tok{2});   % middle of y
        zc = str2double(tok{3});   % middle slice

        % Clamped start positions (1-based)
        x0   = clampStart(xc, W, imgW);
        y0   = clampStart(yc, H, imgH);
        z_lo = clampStart(zc, D, nSlices);
        z_hi = z_lo + min(D, nSlices) - 1;

        fprintf('\n[%d/%d] %s: x0=%d y0=%d z=%d-%d\n', ...
            k, numel(roiNames), roi, x0, y0, z_lo, z_hi);

        outdir = fullfile(outroot, roi);
        cropTiffs(tiffdir, x0, y0, z_lo, z_hi, W, H, outdir);
    end
end

function s = clampStart(center, len, total)
    % Start index of a window of length len centred on `center`,
    % shifted to stay inside [1, total].
    len = min(len, total);
    s = round(center - len/2);
    s = max(s, 1);                 % ran into the start -> take first `len`
    s = min(s, total - len + 1);   % ran into the end   -> take last `len`
end

function cropTiffs(tiffdir, x, y, z_lo, z_hi, w, h, outdir)
    files = [dir(fullfile(tiffdir, '*.tif')); dir(fullfile(tiffdir, '*.tiff'))];
    if isempty(files)
        error('No TIFF files found in directory: %s', tiffdir);
    end
    [~, idx] = sort({files.name});
    names = {files(idx).name};          % plain cell array: cheaper to send to workers

    if ~exist(outdir, 'dir')
        mkdir(outdir);                  % must exist before the parfor
    end

    slices = z_lo:z_hi;
    nOut = numel(slices);

    parfor j = 1:nOut
        name = names{slices(j)};
        I  = imread(fullfile(tiffdir, name));
        Ic = I(y:y+h-1, x:x+w-1, :);
        imwrite(Ic, fullfile(outdir, name));
    end

    fprintf('Finished cropping %d TIFFs into %s\n', nOut, outdir);
end