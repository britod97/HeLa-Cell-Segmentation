function mask = segmentManual_2(inputImage, outputFile)
%SEGMENTMITOCHONDRIAMANUAL Interactive hand-segmentation of mitochondria.
%
%   MASK = SEGMENTMITOCHONDRIAMANUAL(INPUTIMAGE, OUTPUTFILE)
%
%   Loads a 2D image or 3D image volume, lets you draw one or more
%   freehand/assisted ROIs around each mitochondrion on every slice using
%   DRAWASSISTED, and returns (and saves) the union of all drawn ROIs as
%   a single logical mask.
%
%   INPUTIMAGE
%       One of:
%         - Path to a single 2D image file (png/tif/jpg/...)
%         - Path to a multi-page TIFF stack (loaded as a 3D volume)
%         - Path to a folder containing a sequence of 2D image slices
%           (files are sorted alphanumerically and stacked into a volume)
%         - A 2D or 3D numeric array already in the workspace
%
%   OUTPUTFILE (optional)
%       Where to save the result. Saved as a .mat file containing:
%           mask     - logical array, same size as INPUTIMAGE, true
%                      wherever ANY mitochondrion ROI was drawn
%           roiData  - cell array (1 per slice) of cell arrays of ROI
%                      polygon vertices, kept for provenance / re-editing
%       If OUTPUTFILE ends in '.m', the extension is changed to '.mat'
%       automatically (a .m file cannot store array data). If omitted,
%       defaults to 'mitochondria_mask.mat' in the current folder.
%
%   MASK
%       The logical union mask, also returned in the workspace.
%
%   WORKFLOW / CONTROLS
%       For each slice:
%         1. The slice is shown, with any previously accepted ROIs on
%            that slice outlined in yellow.
%         2. Draw a mitochondrion boundary with DRAWASSISTED (click to
%            place vertices / drag to trace, double-click or right-click
%            to close the region).
%         3. A dialog asks whether to add another mitochondrion on this
%            slice, move to the next slice, or finish and save.
%       You can also draw zero ROIs on a slice and immediately move on.
%
%   EXAMPLES
%       % Single 2D image
%       mask = segmentMitochondriaManual('slice120.tif', 'mito_mask.mat');
%
%       % Multi-page TIFF volume
%       mask = segmentMitochondriaManual('EM_stack.tif', 'mito_mask.mat');
%
%       % Folder of slices already cropped (matches CIL50061 workflow)
%       mask = segmentMitochondriaManual('D:\...\Cropped_Tiffs_2', ...
%                                         'mito_mask.mat');
%
%       % Array already in memory
%       mask = segmentMitochondriaManual(myVolume, 'mito_mask.mat');

    if nargin < 2 || isempty(outputFile)
        outputFile = 'mitochondria_mask.mat';
    end
    outputFile = resolveOutputFile(outputFile);

    vol = loadAsVolume(inputImage);
    nSlices = size(vol, 3);
    imgH = size(vol, 1);
    imgW = size(vol, 2);

    mask = false(size(vol));
    roiData = cell(1, nSlices);

    fig = figure('Name', 'Mitochondria Segmentation', ...
                  'NumberTitle', 'off', ...
                  'CloseRequestFcn', @(~,~) finishEarly());
    ax = axes('Parent', fig, 'Units', 'normalized', ...
              'Position', [0.05 0.05 0.62 0.85]);
    % axPrev = axes('Parent', fig, 'Units', 'normalized', ...
    %           'Position', [0.70 0.52 0.28 0.36]);
    % axNext = axes('Parent', fig, 'Units', 'normalized', ...
    %           'Position', [0.70 0.08 0.28 0.36]);

    fig.WindowScrollWheelFcn = @scrollZoomCallback;
    uicontrol('Parent', fig, 'Style', 'pushbutton', ...
        'String', 'Reset View', 'Units', 'normalized', ...
        'Position', [0.01 0.93 0.14 0.05], ...
        'Callback', @resetViewCallback);
    uicontrol('Parent', fig, 'Style', 'text', ...
        'String', 'Scroll to zoom (centered on cursor)  |  magenta = previous slice mask', ...
        'Units', 'normalized', 'Position', [0.16 0.93 0.53 0.05], ...
        'HorizontalAlignment', 'left', 'BackgroundColor', fig.Color);

    stopAll = false;

    for k = 1:nSlices
        if stopAll
            break
        end
        sliceImg = vol(:, :, k);
        sliceRois = {};   % polygon vertex arrays accepted on this slice

        if k > 1
            prevMask = mask(:, :, k - 1);
        else
            prevMask = [];
        end
        % updateNeighborPreviews(axPrev, axNext, vol, k, nSlices);

        keepGoingOnSlice = true;
        while keepGoingOnSlice && ~stopAll
            showSlice(ax, sliceImg, k, nSlices, sliceRois, prevMask);

            roi = drawassisted(ax);

            if ~isvalid(roi) || isempty(roi.Position)
                % User pressed Esc / cancelled without drawing anything
                choice = askNext(true);
            else
                % m = createMask(roi, size(sliceImg, 1), size(sliceImg, 2));
                % mask(:, :, k) = mask(:, :, k) | m;
                % sliceRois{end+1} = roi.Position; %#ok<AGROW>
                % delete(roi);
                % showSlice(ax, sliceImg, k, nSlices, sliceRois, prevMask);
                % choice = askNext(false);
                m = createMask(roi, size(sliceImg, 1), size(sliceImg, 2));

                keepROI = questdlg('Keep this ROI?', ...
                                   'ROI confirmation', ...
                                   'Keep', 'Delete', 'Keep');
                
                if strcmp(keepROI, 'Keep')
                    mask(:, :, k) = mask(:, :, k) | m;
                    sliceRois{end+1} = roi.Position; %#ok<AGROW>
                end
                
                delete(roi);
                
                showSlice(ax, sliceImg, k, nSlices, sliceRois, prevMask);
                
                if strcmp(keepROI, 'Delete')
                    choice = 'another';   % immediately draw again
                else
                    choice = askNext(false);
                end
            end

            switch choice
                case 'another'
                    % loop again on same slice
                case 'nextSlice'
                    keepGoingOnSlice = false;
                case 'finish'
                    keepGoingOnSlice = false;
                    stopAll = true;
            end
        end

        roiData{k} = sliceRois;
    end

    if isvalid(fig)
        delete(fig);
    end

    if ndims(vol) == 2 %#ok<ISMAT>
        mask = mask(:, :, 1);
        roiData = roiData(1);
    end

    save(outputFile, 'mask', 'roiData');
    fprintf('Saved mask (%s) and ROI data to: %s\n', ...
        mat2str(size(mask)), outputFile);

    function finishEarly()
        stopAll = true;
        if isvalid(fig)
            delete(fig);
        end
    end

    function scrollZoomCallback(~, evt)
        if ~isvalid(ax)
            return
        end
        cp = ax.CurrentPoint;
        xCenter = cp(1, 1);
        yCenter = cp(1, 2);

        zoomStep = 1.2;
        if evt.VerticalScrollCount > 0
            f = zoomStep;       % scroll down/back -> zoom out
        else
            f = 1 / zoomStep;   % scroll up/forward  -> zoom in
        end

        xl = xlim(ax);
        yl = ylim(ax);
        newXL = xCenter + ([xl(1), xl(2)] - xCenter) * f;
        newYL = yCenter + ([yl(1), yl(2)] - yCenter) * f;

        % Clamp to image extent so you can't scroll away from the image
        newXL(1) = max(newXL(1), 0.5);
        newXL(2) = min(newXL(2), imgW + 0.5);
        newYL(1) = max(newYL(1), 0.5);
        newYL(2) = min(newYL(2), imgH + 0.5);

        % Ignore if it would zoom in past a sane minimum extent
        if diff(newXL) >= 5 && diff(newYL) >= 5
            xlim(ax, newXL);
            ylim(ax, newYL);
        end
    end

    function resetViewCallback(~, ~)
        if isvalid(ax)
            xlim(ax, [0.5, imgW + 0.5]);
            ylim(ax, [0.5, imgH + 0.5]);
        end
    end
end

% ------------------------------------------------------------------
function outFile = resolveOutputFile(outFile)
    [p, n, ext] = fileparts(outFile);
    if strcmpi(ext, '.m')
        warning(['A .m file cannot store array data. Saving to a ' ...
                  '.mat file with the same name instead.']);
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

    if ~(ischar(inputImage) || isstring(inputImage))
        error('INPUTIMAGE must be a file path, folder path, or numeric array.');
    end
    inputImage = char(inputImage);

    if isfolder(inputImage)
        files = dir(fullfile(inputImage, '*.*'));
        files = files(~[files.isdir]);
        exts = lower(string({files.name}));
        isImg = endsWith(exts, [".tif", ".tiff", ".png", ".jpg", ".jpeg"]);
        files = files(isImg);
        if isempty(files)
            error('No image files found in folder: %s', inputImage);
        end
        [~, order] = sort({files.name});
        files = files(order);

        first = imread(fullfile(inputImage, files(1).name));
        vol = zeros([size(first, 1), size(first, 2), numel(files)], 'like', first);
        vol(:, :, 1) = first;
        for i = 2:numel(files)
            vol(:, :, i) = imread(fullfile(inputImage, files(i).name));
        end
        return
    end

    if ~isfile(inputImage)
        error('File not found: %s', inputImage);
    end

    info = imfinfo(inputImage);
    if numel(info) > 1
        first = imread(inputImage, 1);
        vol = zeros([size(first, 1), size(first, 2), numel(info)], 'like', first);
        vol(:, :, 1) = first;
        for i = 2:numel(info)
            vol(:, :, i) = imread(inputImage, i);
        end
    else
        vol = imread(inputImage);
    end
end

% ------------------------------------------------------------------
function showSlice(ax, sliceImg, k, nSlices, sliceRois, prevMask)
    hadView = ~isempty(ax.Children);
    if hadView
        xl = xlim(ax);
        yl = ylim(ax);
    end

    cla(ax);
    imshow(sliceImg, [], 'Parent', ax);

    if hadView
        xlim(ax, xl);
        ylim(ax, yl);
    end

    hold(ax, 'on');

    hasGhost = ~isempty(prevMask) && any(prevMask(:));
    % if hasGhost
    %     contour(ax, double(prevMask), [0.5, 0.5], ...
    %         'LineColor', [1 0 1], 'LineWidth', 1);
    % end
    if hasGhost
        h = imshow(cat(3, ones(size(prevMask)), zeros(size(prevMask)), ones(size(prevMask))), ...
                   'Parent', ax);     % Magenta RGB image
        h.AlphaData = 0.15 * prevMask;   % 15% opacity only where mask exists
    end

    for i = 1:numel(sliceRois)
        pos = sliceRois{i};
        plot(ax, pos([1:end, 1], 1), pos([1:end, 1], 2), ...
            'y-', 'LineWidth', 1.5);
    end
    hold(ax, 'off');

    if hasGhost
        ghostNote = '  |  magenta = previous slice mask';
    else
        ghostNote = '';
    end
    title(ax, sprintf('Slice %d of %d  |  %d mitochondrion ROI(s) accepted%s', ...
        k, nSlices, numel(sliceRois), ghostNote), 'Interpreter', 'none');
end

% ------------------------------------------------------------------
% function updateNeighborPreviews(axPrev, axNext, vol, k, nSlices)
%     if k > 1
%         imshow(vol(:, :, k - 1), [], 'Parent', axPrev);
%         title(axPrev, sprintf('Slice %d (previous)', k - 1), 'Interpreter', 'none');
%     else
%         cla(axPrev);
%         title(axPrev, 'No previous slice', 'Interpreter', 'none');
%     end
% 
%     if k < nSlices
%         imshow(vol(:, :, k + 1), [], 'Parent', axNext);
%         title(axNext, sprintf('Slice %d (next)', k + 1), 'Interpreter', 'none');
%     else
%         cla(axNext);
%         title(axNext, 'No next slice', 'Interpreter', 'none');
%     end
% end

% ------------------------------------------------------------------
function choice = askNext(emptyDraw)
    if emptyDraw
        msg = 'No ROI drawn. What next?';
        opts = {'Draw again', 'Next slice', 'Finish & save'};
    else
        msg = 'Mitochondrion added. What next?';
        opts = {'Add another', 'Next slice', 'Finish & save'};
    end
    answer = questdlg(msg, 'Continue segmentation', ...
        opts{1}, opts{2}, opts{3}, opts{1});

    switch answer
        case opts{1}
            choice = 'another';
        case opts{2}
            choice = 'nextSlice';
        otherwise
            choice = 'finish';
    end
end
