function mask = segmentManualOctave(inputImage, outputFile)
%SEGMENTMANUAL Interactive hand-segmentation of regions (Octave + MATLAB).
%
%   MASK = SEGMENTMANUALOCTAVE(INPUTIMAGE, OUTPUTFILE)
%
%   INPUTIMAGE: path to a 2D image, a multi-page TIFF, a folder of 2D
%               slices (sorted by name), or a 2D/3D numeric array.
%   OUTPUTFILE: .mat file to save (default 'membrane_mask.mat'),
%               containing  mask (logical, same size as the input) and
%               roiData (1 cell per slice, each a cell of Nx2 [x y] polygons).
%
%   CONTROLS (click inside the figure window first so it has focus)
%       Left click ......... add polygon vertex
%       Right click / Enter  close polygon and accept it as an ROI
%       u or Backspace ..... undo last vertex (or last ROI if no vertices)
%       n .................. next slice (finishes after the last slice)
%       p .................. previous slice
%       + / - .............. zoom in / out, centred on the cursor
%       r .................. reset view
%       f .................. finish and save now
%
%   Requires a graphics toolkit that supports ginput (in Octave:
%   graphics_toolkit('qt') is the usual choice).

    if nargin < 2 || isempty(outputFile)
        outputFile = 'membrane_mask.mat';
    end
    outputFile = resolveOutputFile(outputFile);

    vol = loadAsVolume(inputImage);
    [imgH, imgW, nSlices] = size(vol);

    mask = false(imgH, imgW, nSlices);
    roiData = cell(1, nSlices);
    for i = 1:nSlices
        roiData{i} = {};
    end

    fig = figure('Name', 'Segmentation', 'NumberTitle', 'off');
    ax = axes('Parent', fig);

    k = 1;
    stopAll = false;

    while k <= nSlices && ~stopAll
        hPending = showSlice(ax, vol(:, :, k), k, nSlices, roiData{k});
        verts = zeros(0, 2);
        sliceDone = false;

        while ~sliceDone && ~stopAll
            % ---- get one click / keypress (fails if window is closed) ----
            try
                [x, y, b] = ginput(1);
            catch
                stopAll = true;
                break
            end
            if isempty(b)            % Enter in MATLAB returns empty
                b = 13; x = NaN; y = NaN;
            end

            if b == 1                                   % add vertex
                verts(end+1, :) = [x, y];
                set(hPending, 'XData', verts(:, 1), 'YData', verts(:, 2));

            elseif b == 3 || b == 13 || b == 10         % close polygon
                if size(verts, 1) >= 3
                    roiData{k}{end+1} = verts;
                    mask(:, :, k) = roisToMask(roiData{k}, imgH, imgW);
                end
                verts = zeros(0, 2);
                hPending = showSlice(ax, vol(:, :, k), k, nSlices, roiData{k});

            elseif b == 8 || b == 127 || b == double('u')   % undo
                if ~isempty(verts)
                    verts(end, :) = [];
                    set(hPending, 'XData', verts(:, 1), 'YData', verts(:, 2));
                elseif ~isempty(roiData{k})
                    roiData{k}(end) = [];
                    mask(:, :, k) = roisToMask(roiData{k}, imgH, imgW);
                    hPending = showSlice(ax, vol(:, :, k), k, nSlices, roiData{k});
                end

            elseif b == double('n')                     % next slice
                k = k + 1;
                sliceDone = true;

            elseif b == double('p')                     % previous slice
                k = max(k - 1, 1);
                sliceDone = true;

            elseif b == double('f')                     % finish
                stopAll = true;

            elseif b == double('+') || b == double('=') % zoom in
                zoomAt(ax, x, y, 1/1.3, imgW, imgH);

            elseif b == double('-')                     % zoom out
                zoomAt(ax, x, y, 1.3, imgW, imgH);

            elseif b == double('r')                     % reset view
                set(ax, 'XLim', [0.5, imgW + 0.5], 'YLim', [0.5, imgH + 0.5]);
            end
        end
    end

    if ishghandle(fig)
        delete(fig);
    end

    save('-v7', outputFile, 'mask', 'roiData');
    fprintf('Saved mask (%s) and ROI data to: %s\n', ...
        mat2str(size(mask)), outputFile);
end

% ------------------------------------------------------------------
function outFile = resolveOutputFile(outFile)
    [p, n, ext] = fileparts(outFile);
    if strcmpi(ext, '.m')
        warning('A .m file cannot store array data. Using .mat instead.');
        ext = '.mat';
    elseif isempty(ext)
        ext = '.mat';
    end
    outFile = fullfile(p, [n, ext]);
end

% ------------------------------------------------------------------
function vol = loadAsVolume(inputImage)
    if isnumeric(inputImage) || islogical(inputImage)
        vol = inputImage;
        return
    end
    if ~ischar(inputImage)
        error('INPUTIMAGE must be a file path, folder path, or numeric array.');
    end

    if exist(inputImage, 'dir') == 7
        d = dir(inputImage);
        d = d(~[d.isdir]);
        keep = false(1, numel(d));
        for i = 1:numel(d)
            keep(i) = ~isempty(regexpi(d(i).name, '\.(tif|tiff|png|jpg|jpeg)$', 'once'));
        end
        d = d(keep);
        if isempty(d)
            error('No image files found in folder: %s', inputImage);
        end
        [~, order] = sort({d.name});
        d = d(order);

        first = imread(fullfile(inputImage, d(1).name));
        vol = zeros([size(first, 1), size(first, 2), numel(d)], class(first));
        vol(:, :, 1) = first;
        for i = 2:numel(d)
            vol(:, :, i) = imread(fullfile(inputImage, d(i).name));
        end
        return
    end

    if exist(inputImage, 'file') ~= 2
        error('File not found: %s', inputImage);
    end

    info = imfinfo(inputImage);
    if numel(info) > 1
        first = imread(inputImage, 1);
        vol = zeros([size(first, 1), size(first, 2), numel(info)], class(first));
        vol(:, :, 1) = first;
        for i = 2:numel(info)
            vol(:, :, i) = imread(inputImage, i);
        end
    else
        vol = imread(inputImage);
    end
end

% ------------------------------------------------------------------
function hPending = showSlice(ax, sliceImg, k, nSlices, rois)
    hadView = ~isempty(get(ax, 'Children'));
    if hadView
        xl = get(ax, 'XLim');
        yl = get(ax, 'YLim');
    end

    cla(ax);
    imagesc(ax, double(sliceImg));
    colormap(ax, gray(256));
    axis(ax, 'image');

    if hadView
        set(ax, 'XLim', xl, 'YLim', yl);
    end

    hold(ax, 'on');
    for i = 1:numel(rois)
        pos = rois{i};
        plot(ax, pos([1:end, 1], 1), pos([1:end, 1], 2), 'y-', 'LineWidth', 1.5);
    end
    hPending = plot(ax, NaN, NaN, 'c-o', 'LineWidth', 1, 'MarkerSize', 4);
    hold(ax, 'off');

    title(ax, sprintf(['Slice %d/%d | %d ROI(s) | L-click: vertex, R-click/Enter: close, ' ...
        'u: undo, n/p: next/prev, +/-: zoom, r: reset, f: finish'], ...
        k, nSlices, numel(rois)), 'Interpreter', 'none');
    drawnow;
end

% ------------------------------------------------------------------
function zoomAt(ax, xc, yc, f, imgW, imgH)
    xl = get(ax, 'XLim');
    yl = get(ax, 'YLim');
    newXL = xc + (xl - xc) * f;
    newYL = yc + (yl - yc) * f;
    newXL = [max(newXL(1), 0.5), min(newXL(2), imgW + 0.5)];
    newYL = [max(newYL(1), 0.5), min(newYL(2), imgH + 0.5)];
    if diff(newXL) >= 5 && diff(newYL) >= 5
        set(ax, 'XLim', newXL, 'YLim', newYL);
    end
end

% ------------------------------------------------------------------
function m = roisToMask(rois, imgH, imgW)
    % Union of polygon masks; only tests pixels inside each bounding box.
    m = false(imgH, imgW);
    for i = 1:numel(rois)
        pos = rois{i};
        c1 = max(1, floor(min(pos(:, 1))));  c2 = min(imgW, ceil(max(pos(:, 1))));
        r1 = max(1, floor(min(pos(:, 2))));  r2 = min(imgH, ceil(max(pos(:, 2))));
        if c2 < c1 || r2 < r1
            continue
        end
        [X, Y] = meshgrid(c1:c2, r1:r2);
        m(r1:r2, c1:c2) = m(r1:r2, c1:c2) | inpolygon(X, Y, pos(:, 1), pos(:, 2));
    end
end
