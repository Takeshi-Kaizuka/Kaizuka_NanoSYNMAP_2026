% This script generates tile-wise “heatmap-style” montage images (32-bit single TIFFs)
% using (i) per-tile summary metrics stored in a CSV file and (ii) tile coordinates
% stored in an acquisition metadata text file (*.txt).
% Each tile region in the montage is filled with a constant value for each metric.
%
% Inputs:
% - config/fret_proximity_ratio.json
%     * sample_id: used to locate the input image folder
%     * result_tag: used to search for the summary CSV filename
%     * pixel_size: used to convert metadata coordinates (assumed in µm) to pixels
%     * down_sampling_factor: down-sampling factor for montage generation
% - data/raw_images/{sample_id}/acceptor/
%     * *.txt: metadata file containing tile coordinates (read by read_coordi_metadata)
%     * *.tif: tile images (image content is not used; only height/width are read for sizing/validation)
% - results/summary/
%     * a CSV file whose name contains result_tag (loaded with readtable)
%
% Processing (metrics):
% - Puncta_Density = (Puncta_Number × 100) / 3292.343
% - area_mean is scaled by × 0.004624
% - FRET_efficiency_mean is scaled by × 100 (percentage)
% - other metrics are used as-is from the CSV
% - tiles missing from the CSV are filled with 0
% - standard-deviation (“_std”) metrics are not handled in this version
%
% Outputs:
% - results/tile_heatmap/{sample_id}/
%   32-bit float TIFFs saved as "Montage_{folder_label}_{metric}.tif" for each metric.

%% 1. Initial setup (relative to the project root)
script_dir   = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);

config_path = fullfile(project_root, 'config', 'fret_proximity_ratio.json');
cfg = jsondecode(fileread(config_path));

sample_id  = cfg.sample_id;
result_tag = cfg.result_tag;

if ~isfield(cfg, 'pixel_size') || ~isfield(cfg, 'down_sampling_factor')
    error('Config file must contain "pixel_size" and "down_sampling_factor": %s', config_path);
end

pixel_size = cfg.pixel_size;
down_sampling_factor = cfg.down_sampling_factor;

% Input folder: data/raw_images/{sample_id}/acceptor/
acceptor_dir = fullfile(project_root, 'data', 'raw_images', sample_id, 'acceptor');
if ~isfolder(acceptor_dir)
    error('Input folder not found: %s', acceptor_dir);
end

select_mouse_path = {acceptor_dir};

% Output folder: results/tile_heatmap/{sample_id}/
out_dir = fullfile(project_root, 'results', 'tile_heatmap', sample_id);
if ~isfolder(out_dir)
    mkdir(out_dir);
end

%% 2. Load the CSV file (search in results/summary using result_tag)
summary_dir = fullfile(project_root, 'results', 'summary');
if ~isfolder(summary_dir)
    error('Summary folder not found: %s', summary_dir);
end

csv_files = dir(fullfile(summary_dir, '*.csv'));
target_csv = '';

for k = 1:length(csv_files)
    if contains(csv_files(k).name, result_tag)
        target_csv = csv_files(k).name;
        break;
    end
end

if isempty(target_csv)
    error('No CSV containing result_tag "%s" was found in: %s', result_tag, summary_dir);
end

csv_file = fullfile(summary_dir, target_csv);
summary_data = readtable(csv_file);


%% 3. Generate montage images for each target folder
metrics = {'Puncta_Density', 'area_mean', 'inverse_aspect_mean', 'Donor_mean_intensity_mean', 'Acceptor_mean_intensity_mean', 'FRET_mean_intensity_mean', 'FRET_efficiency_mean', 'CC_Acceptor_FRETeff'};

tic
for i = 1:length(select_mouse_path)
    % Retrieve metadata (*.txt) and images (*.tif)
    dir_meta = dir(fullfile(select_mouse_path{i}, '*.txt'));
    dir_tif = dir(fullfile(select_mouse_path{i}, '*.tif'));
    
    if ~isempty(dir_meta)
        % Read tile coordinates from the metadata
        coordi_xy = read_coordi_metadata(fullfile(select_mouse_path{i}, dir_meta(1).name), length(dir_tif), pixel_size);
        
        % Get the image size from the first TIFF image
        iminfo = imfinfo(fullfile(select_mouse_path{i}, dir_tif(1).name));
        im_height = iminfo(1).Height;
        im_width = iminfo(1).Width;
        im_height_downsample = ceil(im_height / down_sampling_factor) + 1;
        im_width_downsample = ceil(im_width / down_sampling_factor) + 1;
        
        % Compute downsampled coordinates and determine the overall montage size
        coordi_xy_downsample = round(coordi_xy / down_sampling_factor);
        Montage_dim = max(coordi_xy_downsample);
        
        % Initialize montage images for each metric as 32-bit (single)
        Montage = struct();
        for m = 1:length(metrics)
            Montage.(metrics{m}) = nan(Montage_dim(2) + im_height_downsample + 1, Montage_dim(1) + im_width_downsample + 1, 'single');
        end
        
        % Process each tile (find the row where CSV 'File' matches i_img)
        for i_img = 1:length(dir_tif)
            temp_xy_downsample = coordi_xy_downsample(i_img, :);
            row_index = find(summary_data.File == i_img, 1);
            if ~isempty(row_index)
                % Access each metric value by its header name and perform the required calculations
                puncta_number = summary_data.Puncta_Number(row_index);
                puncta_density = (puncta_number * 100) / 3292.343;          
                area_mean_val = summary_data.area_mean(row_index) * 0.004624;                
                inverse_aspect_mean_val = summary_data.inverse_aspect_mean(row_index);                
                donor_mean_intensity_mean_val = summary_data.Donor_mean_intensity_mean(row_index);               
                acceptor_mean_intensity_mean_val = summary_data.Acceptor_mean_intensity_mean(row_index);               
                fret_mean_intensity_mean_val = summary_data.FRET_mean_intensity_mean(row_index);               
                fret_efficiency_mean_val = summary_data.FRET_efficiency_mean(row_index) * 100;                
                cc_acceptor_freteff_val = summary_data.CC_Acceptor_FRETeff(row_index);
            else
                % If no corresponding data exist, set all values to 0
                puncta_density = 0;
                area_mean_val = 0;
                inverse_aspect_mean_val = 0;
                donor_mean_intensity_mean_val = 0;
                acceptor_mean_intensity_mean_val = 0;
                fret_mean_intensity_mean_val = 0;
                fret_efficiency_mean_val = 0;
                cc_acceptor_freteff_val = 0;
            end
            
            % Compute the tile placement position
            y_start = temp_xy_downsample(2) + 1;
            y_end = temp_xy_downsample(2) + im_height_downsample;
            x_start = temp_xy_downsample(1) + 1;
            x_end = temp_xy_downsample(1) + im_width_downsample;
            
            % Assign the tile values into each montage image
            Montage.Puncta_Density(y_start:y_end, x_start:x_end) = puncta_density;
            Montage.area_mean(y_start:y_end, x_start:x_end) = area_mean_val;
            Montage.inverse_aspect_mean(y_start:y_end, x_start:x_end) = inverse_aspect_mean_val;
            Montage.Donor_mean_intensity_mean(y_start:y_end, x_start:x_end) = donor_mean_intensity_mean_val;
            Montage.Acceptor_mean_intensity_mean(y_start:y_end, x_start:x_end) = acceptor_mean_intensity_mean_val;
            Montage.FRET_mean_intensity_mean(y_start:y_end, x_start:x_end) = fret_mean_intensity_mean_val;
            Montage.FRET_efficiency_mean(y_start:y_end, x_start:x_end) = fret_efficiency_mean_val;
            Montage.CC_Acceptor_FRETeff(y_start:y_end, x_start:x_end) = cc_acceptor_freteff_val;
        end
        
        % Set NaN regions to the background value (0)
        for m = 1:length(metrics)
            Montage.(metrics{m})(isnan(Montage.(metrics{m}))) = 0;
        end
        
        % Create the output filename (based on the folder name)
        temp_foldername = strsplit(select_mouse_path{i}, filesep);
        folder_label = temp_foldername{end};
        
        % Write each montage image using the Tiff class (save as single precision)
		for m = 1:length(metrics)
			output_filename = fullfile(out_dir, sprintf('Montage_%s_%s_%s.tif', folder_label, metrics{m}, sample_id));
			writeSingleTiff(Montage.(metrics{m}), output_filename);
		end

        t = toc;
        disp(sprintf('Montage images of %s generated in %.2f seconds.', folder_label, t));
    else
        disp('Metadata cannot be found.');
    end
end

%% Helper function
function coordi_xy = read_coordi_metadata(file_path, num_img, pixel_size)
    fid = fopen(file_path, 'r', 'n', 'Unicode');
    TextAsCells = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    ind = find(~cellfun(@isempty, strfind(TextAsCells{1}, 'Point Name	X Pos[')));
    if isempty(ind)
        error('Cannot find coordinate info in the metadata file.');
        coordi_xy = [];
    else
        coordi_um = zeros(num_img, 2);
        for i_line = ind + 1 : ind + num_img
            temp_text_split = strsplit(TextAsCells{1}{i_line}, '\t');
            temp_locatehash = strfind(temp_text_split{1}, '#');
            temp_tilenum = str2double(temp_text_split{1}(temp_locatehash+1:end));
            temp_x_um = str2double(temp_text_split{2});
            temp_y_um = str2double(temp_text_split{3});
            coordi_um(temp_tilenum, 1) = temp_x_um;
            coordi_um(temp_tilenum, 2) = temp_y_um;
        end
        min_coordi = min(coordi_um);
        coordi_um(:,1) = coordi_um(:,1) - min_coordi(1);
        coordi_um(:,2) = coordi_um(:,2) - min_coordi(2);
        coordi_xy = (coordi_um / pixel_size) + 1;
    end
end

function writeSingleTiff(I, filename)
    t = Tiff(filename, 'w');
    tagstruct.ImageLength = size(I, 1);
    tagstruct.ImageWidth = size(I, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    t.setTag(tagstruct);
    t.write(I);
    t.close();
end
