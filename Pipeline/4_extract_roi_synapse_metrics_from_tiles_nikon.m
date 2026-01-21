% This script extracts synapse metrics (tile-level “Summary_metrics_*.csv”) that fall within ROIs,
% and exports ROI-wise aggregated CSV files.
%
% [Inputs / referenced files]
% - Config : config/fret_proximity_ratio.json
%           The script reads sample_id, down_sampling_factor, and pixel_size from this file.
% - Images : TIFF tiles (*.tif) and metadata (*.txt) in data/raw_images/{sample_id}/acceptor/
% - CSV    : Tile-wise synapse metrics, Summary_metrics_*.csv, in data/csv/{sample_id}/
% - ROIs   : ROI files (*.roi) or an ROI.zip in results/roi_data/roi/
%           If {sample_id}.zip exists, it is used preferentially.
%
% [Outputs]
% - results/roi_data/{sample_id}/
%     - ROI-wise aggregated results: [ROI_name].csv
%     - List of ROIs skipped due to read errors: SkippedFiles.txt
%
% [Workflow] (for all ROIs in results/roi_data/roi/)
% 1) Build input/output paths using sample_id from the config and prepare the output folder
%    (results/roi_data/{sample_id}/).
% 2) From path1 (TIFF tiles + metadata .txt), read stage coordinates for each tile,
%    convert them to pixel coordinates using pixel_size, then convert to the downsampled coordinate system
%    using down_sampling_factor to compute each tile’s rectangular region (tile position).
% 3) Load ROIs from path3 (results/roi_data/roi/; either *.roi or *.roi inside a zip) and transform ROI
%    coordinates back to the analysis coordinate system (the script assumes the inverse transform described
%    in comments, e.g., “90° clockwise rotation + left-right flip” via coordinate swapping, etc.).
% 4) For each ROI, test overlaps with all tile rectangles and extract rows from tile CSVs:
%    - If the full tile is inside the ROI: include all rows from that tile CSV.
%    - If partially overlapping: convert each row’s centroid to the downsampled montage coordinates and
%      include only rows whose centroids fall within the ROI.
% 5) Replace the first column of each included data row with the tile number parsed from the CSV filename.
% 6) Export the aggregated table for each ROI as [ROI_name].csv into results/roi_data/{sample_id}/.
% 7) Skip ROIs that cause read errors (e.g., invalid coordinate counts), and report/save the skipped list
%    as SkippedFiles.txt.

%% 1. Set paths and parameters
thisFile   = mfilename('fullpath');
scriptsDir = fileparts(thisFile);
projectRoot = fileparts(scriptsDir);

configPath = fullfile(projectRoot, 'config', 'fret_proximity_ratio.json');
if ~isfile(configPath)
    error('Config file not found: %s', configPath);
end
cfg = jsondecode(fileread(configPath));
sample_id = string(cfg.sample_id);

% path1: TIFF images and metadata (.txt)
path1 = fullfile(projectRoot, 'data', 'raw_images', sample_id, 'acceptor');

% path2: Tile-level synapse metrics CSV files
path2 = fullfile(projectRoot, 'data', 'csv', sample_id);

% path3: ROI data
roi_input_dir = fullfile(projectRoot, 'results', 'roi_data','roi');
path3 = roi_input_dir;

% Output folder
output_dir = fullfile(projectRoot, 'results', 'roi_data', sample_id);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Load parameters from the config file
if ~isfield(cfg, 'down_sampling_factor')
    error('Config missing field "down_sampling_factor": %s', configPath);
end
if ~isfield(cfg, 'pixel_size')
    error('Config missing field "pixel_size": %s', configPath);
end

down_sampling_factor = double(cfg.down_sampling_factor);
pixel_size = double(cfg.pixel_size);

% Validate parameter values
if ~isfinite(down_sampling_factor) || down_sampling_factor <= 0 || mod(down_sampling_factor, 1) ~= 0
    error('down_sampling_factor must be a positive integer. Got: %g', down_sampling_factor);
end
if ~isfinite(pixel_size) || pixel_size <= 0
    error('pixel_size must be a positive number. Got: %g', pixel_size);
end

%% 2. Compute tile positions from TIFF images/metadata and apply
tifFiles = dir(fullfile(path1, '*.tif'));
numTiles = numel(tifFiles);
if numTiles == 0
    error('No TIFF files were found in path1.');
end

imInfo = imfinfo(fullfile(path1, tifFiles(1).name));
im_height = imInfo(1).Height;
im_width  = imInfo(1).Width;
im_height_down = ceil(im_height / down_sampling_factor) + 1;
im_width_down  = ceil(im_width / down_sampling_factor) + 1;

metaFiles = dir(fullfile(path1, '*.txt'));
if isempty(metaFiles)
    error('No metadata (.txt) files were found in path1.');
end
metaFilePath = fullfile(path1, metaFiles(1).name);

coordi_xy = read_coordi_metadata(metaFilePath, numTiles, pixel_size);
coordi_xy_down = round(coordi_xy / down_sampling_factor);

tileBoxes = zeros(numTiles, 4);
for i = 1:numTiles
    x = coordi_xy_down(i, 1);
    y = coordi_xy_down(i, 2);
    x_min = x + 1;
    y_min = y + 1;
    x_max = x + im_width_down;
    y_max = y + im_height_down;
    tileBoxes(i, :) = [x_min, y_min, x_max, y_max];
end

%% 3. Process ROI files
% Preferentially search for a ZIP file that matches sample_id
zipFiles = dir(fullfile(path3, sample_id + ".zip"));

if ~isempty(zipFiles)
    zipFilePath = fullfile(path3, zipFiles(1).name);
    tempRoiFolder = fullfile(tempdir, 'Extracted_ROI');
    if exist(tempRoiFolder, 'dir')
        rmdir(tempRoiFolder, 's');
    end
    mkdir(tempRoiFolder);
    unzip(zipFilePath, tempRoiFolder);
    roiFolder = tempRoiFolder;
else
    roiFolder = path3;  % If no ZIP is found, read *.roi files directly under path3
end

roiFiles = dir(fullfile(roiFolder, '*.roi'));
if isempty(roiFiles)
    error('No ROI files were found.');
end

skippedFiles = {};

%% 4. Extract synapse data for each ROI file
for r = 1:length(roiFiles)
    roiFilePath = fullfile(roiFolder, roiFiles(r).name);
    
    try
        roi = readImageJROI(roiFilePath);
    catch ME
        if contains(ME.message, 'The number of coordinates read is invalid.')
            warning('Skipping ROI file %s due to an invalid number of coordinates.', roiFiles(r).name);
            skippedFiles{end+1} = roiFiles(r).name;
            continue;
        else
            rethrow(ME);
        end
    end
    
    % Create an ROI polygon
    % Apply the inverse transform (rotate 90° CW and mirror) to restore the original coordinate system
    if isfield(roi, 'subROIs')
        combinedPoly = [];
        for i = 1:length(roi.subROIs)
            subROI = roi.subROIs{i};
            if numel(subROI.x) < 3 || numel(subROI.y) < 3
                warning('Sub-ROI %d does not have enough vertices.', i);
                continue;
            end
            subPoly = polyshape(subROI.y, subROI.x);
            if isempty(combinedPoly)
                combinedPoly = subPoly;
            else
                combinedPoly = union(combinedPoly, subPoly);
            end
        end
        if isempty(combinedPoly)
            warning('No valid sub-ROIs were found in ROI file %s.', roiFiles(r).name);
            continue;
        end
        roiPoly = combinedPoly;
    else
        if numel(roi.x) < 3 || numel(roi.y) < 3
            warning('ROI file %s does not have enough vertices.', roiFiles(r).name);
            continue;
        end
        roiPoly = polyshape(roi.y, roi.x);
    end
    
    %% 5. Determine overlap with tile CSVs and extract synapse data
    csvFiles = dir(fullfile(path2, 'Summary_metrics_*.csv'));
    if isempty(csvFiles)
        error('No CSV files were found in path2.');
    end
    
    combinedData = [];
    
    for k = 1:length(csvFiles)
        csvFileName = csvFiles(k).name;
        parts = strsplit(csvFileName, '_');
        tileNumStr = parts{end};
        tileNumStr = erase(tileNumStr, '.csv');
        tileNum = str2double(tileNumStr);
        if isnan(tileNum)
            warning('Failed to extract a tile number from filename %s.', csvFileName);
            continue;
        end
        
        if tileNum < 1 || tileNum > numTiles
            warning('Tile number %d in file %s is out of range.', csvFileName, tileNum);
            continue;
        end
        tileBox = tileBoxes(tileNum, :);
        rect_x = [tileBox(1), tileBox(3), tileBox(3), tileBox(1)];
        rect_y = [tileBox(2), tileBox(2), tileBox(4), tileBox(4)];
        tilePoly = polyshape(rect_x, rect_y);
        
        interPoly = intersect(roiPoly, tilePoly);
        if area(interPoly) == 0
            continue;
        end
        
        [in, on] = inpolygon(rect_x, rect_y, roiPoly.Vertices(:,1), roiPoly.Vertices(:,2));
        fullyInside = all(in | on);
        
        T = readtable(fullfile(path2, csvFileName));
        if height(T) < 1
            continue;
        end
        
        if fullyInside
            includeRows = true(height(T), 1);
        else
            xOffset = coordi_xy_down(tileNum, 1);
            yOffset = coordi_xy_down(tileNum, 2);
            centroidsX = T{:,3} / down_sampling_factor + xOffset;
            centroidsY = T{:,4} / down_sampling_factor + yOffset;
            [inPoints, onPoints] = inpolygon(centroidsX, centroidsY, roiPoly.Vertices(:,1), roiPoly.Vertices(:,2));
            includeRows = (inPoints | onPoints);
        end
        
        if any(includeRows)
            T_included = T(includeRows, :);
            T_included{:,1} = repmat(tileNum, height(T_included), 1);
            combinedData = [combinedData; T_included];
        end
    end
    
    %% 6. Export integrated data to CSV
    outputHeader = {'TileNumber','area','centroid-X','centroid-Y','major_axis_length',...
        'minor_axis_length','Donor_mean_intensity','Acceptor_mean_intensity','FRET_mean_intensity','FRET_efficiency'};
    
    if isempty(combinedData)
        warning('No synapse data were found within ROI file %s.', roiFiles(r).name);
    else
        combinedData.Properties.VariableNames = outputHeader;
        [~, roiName, ~] = fileparts(roiFilePath);
        outputCSVFile = fullfile(output_dir, [roiName, '.csv']);
        writetable(combinedData, outputCSVFile);
        fprintf('Saved synapse data within ROI %s to %s. \n', roiFiles(r).name, outputCSVFile);
    end
end

%% 7. Report skipped ROI files and delete the temporary folder
skippedFilePath = fullfile(output_dir, 'SkippedFiles.txt');
fid = fopen(skippedFilePath, 'w');
if ~isempty(skippedFiles)
    fprintf('The following ROI files were skipped due to invalid coordinate information: \n');
    fprintf(fid, 'The following ROI files were skipped due to invalid coordinate information: \n');
    for i = 1:length(skippedFiles)
        fprintf('  %s\n', skippedFiles{i});
        fprintf(fid, '  %s\n', skippedFiles{i});
    end
else
    fprintf('No ROI files were skipped. \n');
    fprintf(fid, 'No ROI files were skipped. \n');
end
fclose(fid);

if ~strcmp(roiFolder, path3)
    rmdir(roiFolder, 's');
end

%% --- Subfunctions ---
function coordi_xy = read_coordi_metadata(file_path, num_img, pixel_size)
    fid = fopen(file_path, 'r');
    if fid == -1
        error('Failed to open metadata file %s.', file_path);
    end
    textLines = textscan(fid, '%s', 'Delimiter', '\n');
    textLines = textLines{1};
    fclose(fid);
    
    ind = find(~cellfun(@isempty, strfind(textLines, 'Point Name	X Pos[')));
    if isempty(ind)
        error('No coordinate information was found in the metadata file.');
    end
    coordi_um = zeros(num_img, 2);
    for i = 1:num_img
        line = textLines{ind + i};
        parts = strsplit(line, '\t');
        pos = strfind(parts{1}, '#');
        tileNum = str2double(parts{1}(pos+1:end));
        x_um = str2double(parts{2});
        y_um = str2double(parts{3});
        coordi_um(tileNum, :) = [x_um, y_um];
    end
    min_coord = min(coordi_um, [], 1);
    coordi_um(:,1) = coordi_um(:,1) - min_coord(1);
    coordi_um(:,2) = coordi_um(:,2) - min_coord(2);
    coordi_xy = (coordi_um / pixel_size) + 1;
end

function roi = readImageJROI(filename)
    fid = fopen(filename, 'r', 'ieee-be');
    if fid == -1
        error('Failed to open ROI file %s.', filename);
    end
    
    magic = fread(fid, 4, 'char')';
    if ~strcmp(char(magic), 'Iout')
        fclose(fid);
        error('This is not an ImageJ ROI file.');
    end
    version = fread(fid, 1, 'int16');
    roiType = fread(fid, 1, 'uint8');
    options = fread(fid, 1, 'uint8');
    top    = fread(fid, 1, 'int16');
    left   = fread(fid, 1, 'int16');
    bottom = fread(fid, 1, 'int16');
    right  = fread(fid, 1, 'int16');
    
    % For a non-composite ROI, read nCoordinates from the header
    isComposite = bitand(options, 128) > 0;
    if ~isComposite
        nCoordinates = fread(fid, 1, 'uint16');
        if isempty(nCoordinates) || nCoordinates <= 0
            fclose(fid);
            error('The number of coordinates read is invalid.: %d', nCoordinates);
        end
        if nCoordinates > 10000
            warning('The number of coordinates read exceeds 10,000.: %d. Reading all coordinates.', nCoordinates);
        end
    end
    
    % The header has a fixed length of 64 bytes
    fseek(fid, 64, 'bof');
    
    if ~isComposite
        % If the subpixel option is enabled, skip an additional 4 bytes (two float values)
        if bitand(options, 8) > 0
            fseek(fid, 4, 'cof');
        end
        x_coords = fread(fid, nCoordinates, 'int16');
        y_coords = fread(fid, nCoordinates, 'int16');
        roi.x = double(x_coords) + double(left);
        roi.y = double(y_coords) + double(top);
    else
        % Handle composite ROI
        nSubROIs = fread(fid, 1, 'int16');
        if isempty(nSubROIs) || nSubROIs <= 0
            fclose(fid);
            error('Invalid number of sub-ROIs in the composite ROI.: %d', nSubROIs);
        end
        roi.subROIs = cell(nSubROIs, 1);
        for i = 1:nSubROIs
            nCoord = fread(fid, 1, 'uint16');
            if isempty(nCoord) || nCoord <= 0
                fclose(fid);
                error('Invalid number of vertices in sub-ROI %d.: %d', i, nCoord);
            end
            x_coords = fread(fid, nCoord, 'int16');
            y_coords = fread(fid, nCoord, 'int16');
            subROI.x = double(x_coords) + double(left);
            subROI.y = double(y_coords) + double(top);
            roi.subROIs{i} = subROI;
        end
    end
    
    fclose(fid);
end
