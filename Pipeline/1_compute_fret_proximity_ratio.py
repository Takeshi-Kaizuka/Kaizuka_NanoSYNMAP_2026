#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script matches the TIFF image files (donor, acceptor, and FRET channels) in each folder,
performs background correction, fluorescence leakage correction, FRET efficiency calculation,
and outputs the results to CSV files.
"""

import json
import os
import sys
import glob
from pathlib import Path
import numpy as np
import pandas as pd
from skimage.io import imread
from skimage import filters, measure
from skimage.morphology import remove_small_objects

# =============================
# Specify folder paths (project-root based)

def find_project_root(start_dir: Path) -> Path:
    """
    Find project root by looking for README.md (and optionally environment.yml).
    This makes the script robust even if it is executed from a different working directory.
    """
    for p in [start_dir] + list(start_dir.parents):
        if (p / "README.md").exists():
            return p
    return start_dir

# project-root
project_root = find_project_root(Path(__file__).resolve().parent)

# =============================
# Load config (JSON) from project-root/config/
CONFIG_PATH = project_root / "config" / "fret_proximity_ratio.json"
if not CONFIG_PATH.exists():
    raise FileNotFoundError(f"Config file not found: {CONFIG_PATH}")

with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    cfg = json.load(f)

# Required parameters
sample_id = cfg["sample_id"]
result_tag = cfg["result_tag"]
start_index = int(cfg["start_index"])

SATURATION_VALUE = int(cfg["saturation_value"])
MIN_SIZE = int(cfg["mask_min_size"])
MAX_AREA = int(cfg["mask_max_area"])

DONOR_BLEEDTHROUGH = float(cfg["donor_bleedthrough"])
ACCEPTOR_CROSSTALK = float(cfg["acceptor_crosstalk"])

# =============================
donor_folder    = project_root / "data" / "raw_images" / sample_id / "donor"
acceptor_folder = project_root / "data" / "raw_images" / sample_id / "acceptor"
fret_folder     = project_root / "data" / "raw_images" / sample_id / "fret"
output_folder   = project_root / "data" / "csv" / sample_id
output_folder.mkdir(parents=True, exist_ok=True)
# =============================

def load_image(filepath):
    """Load an image from a file path."""
    return imread(filepath)

def threshold_image_otsu(input_image):
    """Perform Otsu thresholding on the input image and return a binary image.
    Pixels with np.nan are excluded from the threshold calculation and set to False.
    """
    valid_pixels = input_image[np.isfinite(input_image)]
    threshold_value = filters.threshold_otsu(valid_pixels)
    binary = np.zeros_like(input_image, dtype=bool)
    binary[np.isfinite(input_image)] = input_image[np.isfinite(input_image)] > threshold_value
    return binary

def coincidence_analysis_pixels(binary_image1, binary_image2):
    """Calculate the logical AND of two binary images and return the overlapping mask."""
    return binary_image1 & binary_image2

def analyse_labelled_image(labelled_image, original_image):
    """
    Analyze labelled regions from the image.
    Returns a DataFrame containing area, centroid, major axis length,
    minor axis length, and mean intensity.
    """
    props = measure.regionprops_table(
        labelled_image, intensity_image=original_image,
        properties=('area','centroid','major_axis_length','minor_axis_length','mean_intensity')
    )
    return pd.DataFrame.from_dict(props)

def label_image(input_image):
    """Label connected regions in a binary image."""
    return measure.label(input_image)

def remove_saturated_regions(original_image, binary_mask, saturation_value):
    """
    For each connected region in the binary_mask obtained from initial thresholding,
    check the corresponding region in the original_image. If the region contains any
    saturated pixels (value 4095), set all pixel values in that region to np.nan.
    Return the modified image.
    """
    labelled = measure.label(binary_mask)
    new_image = original_image.copy().astype(float)
    for region in measure.regionprops(labelled, intensity_image=original_image):
        region_pixels = original_image[labelled == region.label]
        if np.any(region_pixels == saturation_value):
            new_image[labelled == region.label] = np.nan
    return new_image

# Get list of TIFF files from each folder
donor_files    = sorted(glob.glob(os.path.join(donor_folder, "*.tif")))
acceptor_files = sorted(glob.glob(os.path.join(acceptor_folder, "*.tif")))
fret_files     = sorted(glob.glob(os.path.join(fret_folder, "*.tif")))

# Check that the number of files in each folder is the same
num_donor    = len(donor_files)
num_acceptor = len(acceptor_files)
num_fret     = len(fret_files)

if not (num_donor == num_acceptor == num_fret):
    sys.exit("Error: The number of TIFF files in each folder does not match. Donor: {}, Acceptor: {}, FRET: {}"
             .format(num_donor, num_acceptor, num_fret))
else:
    print("The number of image files in the three channels has been confirmed to match.")

# Check that the specified start index is within range (1-indexed)
if start_index < 1 or start_index > num_donor:
    sys.exit("Error: The specified start index {} is out of range.".format(start_index))

# Process each group of files (from start_index to the last image in the sorted list)
for idx in range(start_index - 1, num_donor):
    donor_file    = donor_files[idx]
    acceptor_file = acceptor_files[idx]
    fret_file     = fret_files[idx]
    
    print("Processing file group {}:".format(idx + 1))
    print("  Donor:    {}".format(os.path.basename(donor_file)))
    print("  Acceptor: {}".format(os.path.basename(acceptor_file)))
    print("  FRET:     {}".format(os.path.basename(fret_file)))
    
    # Load images
    Donor_image    = load_image(donor_file)
    Acceptor_image = load_image(acceptor_file)
    FRET_image     = load_image(fret_file)
    
    # Check for saturated pixels (value 4095) in Donor and Acceptor images
    if np.any(Donor_image == SATURATION_VALUE) or np.any(Acceptor_image == SATURATION_VALUE):
        print("  Saturated pixels detected; applying modified thresholding procedure.")
        # Donor image: initial thresholding, exclude regions with saturated pixels, and re-threshold
        donor_binary_initial = threshold_image_otsu(Donor_image)
        Donor_modified = remove_saturated_regions(Donor_image, donor_binary_initial, SATURATION_VALUE)
        donor_binary_image = threshold_image_otsu(Donor_modified)
        
        # Acceptor image: initial thresholding, exclude regions with saturated pixels, and re-threshold
        acceptor_binary_initial = threshold_image_otsu(Acceptor_image)
        Acceptor_modified = remove_saturated_regions(Acceptor_image, acceptor_binary_initial, SATURATION_VALUE)
        acceptor_binary_image = threshold_image_otsu(Acceptor_modified)
    else:
        print("  No saturated pixels detected; proceeding with standard thresholding.")
        donor_binary_image = threshold_image_otsu(Donor_image)
        acceptor_binary_image = threshold_image_otsu(Acceptor_image)
    
    # Extract overlapping region from the thresholded images
    pixel_overlap_image = coincidence_analysis_pixels(donor_binary_image, acceptor_binary_image)
    
    # Remove small objects (regions smaller than 5 pixels)
    pixel_overlap_image = remove_small_objects(pixel_overlap_image, min_size=MIN_SIZE)
    
    # Remove large regions (exclude regions larger than 5000 pixels)
    labelled_temp = measure.label(pixel_overlap_image)
    props = measure.regionprops(labelled_temp)
    mask_filtered = np.zeros_like(pixel_overlap_image, dtype=bool)
    for region in props:
        if region.area <= MAX_AREA:
            mask_filtered[labelled_temp == region.label] = True
    pixel_overlap_image = mask_filtered
    
    # Background correction using the original images
    donor_bg    = Donor_image.min()
    acceptor_bg = Acceptor_image.min()
    fret_bg     = FRET_image.min()
    
    Donor_corrected    = Donor_image - donor_bg
    Acceptor_corrected = Acceptor_image - acceptor_bg
    FRET_bg_corrected  = FRET_image - fret_bg
    
    # Fluorescence leakage correction
    FRET_corrected = (
     FRET_bg_corrected
     - (DONOR_BLEEDTHROUGH * Donor_corrected)
     - (ACCEPTOR_CROSSTALK * Acceptor_corrected)
    )
   
    # Calculation of Proximity Ratio (FRET efficiency)
    Proximity_image = FRET_corrected / (FRET_corrected + Donor_corrected)
    
    # Labelling and region analysis for each image
    labelled_coincident_image = label_image(pixel_overlap_image)
    FRET_analysis    = analyse_labelled_image(labelled_coincident_image, Proximity_image)
    don_analysis     = analyse_labelled_image(labelled_coincident_image, Donor_corrected)
    acc_analysis     = analyse_labelled_image(labelled_coincident_image, Acceptor_corrected)
    FRETintensity_analysis = analyse_labelled_image(labelled_coincident_image, FRET_corrected)
    
    # Rename columns for centroid coordinates (don_analysis provides centroid-1 as X and centroid-0 as Y)
    don_analysis = don_analysis.rename(columns={"centroid-1": "centroid-X", "centroid-0": "centroid-Y"})
    
    # Compile results and output to CSV
    summary_df = pd.DataFrame({
        "label": np.arange(len(don_analysis)),
        "area": don_analysis["area"],
        "centroid-X": don_analysis["centroid-X"],
        "centroid-Y": don_analysis["centroid-Y"],
        "major_axis_length": don_analysis["major_axis_length"],
        "minor_axis_length": don_analysis["minor_axis_length"],
        "Donor_mean_intensity": don_analysis["mean_intensity"],
        "Acceptor_mean_intensity": acc_analysis["mean_intensity"],
        "FRET_mean_intensity": FRETintensity_analysis["mean_intensity"],
        "FRET_efficiency": FRET_analysis["mean_intensity"]
    })
    out_summary = output_folder / f"Summary_metrics_{result_tag}_{idx + 1}.csv"
    summary_df.to_csv(out_summary, index=False)
    
    print(f"Finished processing group {idx + 1} -> {out_summary.name}")

