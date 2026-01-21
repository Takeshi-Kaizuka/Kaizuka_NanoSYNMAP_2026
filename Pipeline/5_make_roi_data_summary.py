# -*- coding: utf-8 -*-

import os
import glob
import json
from pathlib import Path
import pandas as pd

# ------------------------------------------------------------
# Use the project root as the base path and refer to results/roi_data/{sample_id}/
#  - This script is assumed to be located in scripts/
#  - sample_id is read from "sample_id" in config/fret_proximity_ratio.json
# ------------------------------------------------------------

try:
    project_root = Path(__file__).resolve().parent.parent
except NameError:
    project_root = Path.cwd()

# Load configuration
config_path = project_root / "config" / "fret_proximity_ratio.json"
with open(config_path, "r", encoding="utf-8") as f:
    config = json.load(f)

sample_id = config["sample_id"]

# Input folder (results/roi_data/{sample_id}/)
folder_path = project_root / "results" / "roi_data" / sample_id
if not folder_path.exists():
    raise FileNotFoundError(f"Input folder not found: {folder_path}")

# Collect paths to all CSV files in the folder
csv_files = [p for p in folder_path.glob("*.csv") if p.name != "ROI_Data_Summary.csv"]

# List to store summary results
summary_list = []

# Scale factor for area
scale_factor = 0.068 ** 2

for file in csv_files:
    roi_name = file.stem
    df = pd.read_csv(file)

    # Calculations
    # 1. area: multiply the mean by scale_factor
    area_mean = df['area'].mean() * scale_factor

    # 2. Aspect ratio: mean of (minor_axis_length / major_axis_length) computed row-wise
    ratio = df['minor_axis_length'] / df['major_axis_length']
    aspect_ratio_mean = ratio.mean()
    
    # 3. Mean of Donor_mean_intensity
    donor_mean = df['Donor_mean_intensity'].mean()
    
    # 4. Mean of Acceptor_mean_intensity
    acceptor_mean = df['Acceptor_mean_intensity'].mean()
    
    # 5. Mean of FRET_mean_intensity
    fret_mean = df['FRET_mean_intensity'].mean()
    
    # 6. Mean of FRET_efficiency
    fret_eff_mean = df['FRET_efficiency'].mean() * 100
    
    # Store results as a dictionary
    summary_list.append({
         "ROI": roi_name,
         "area": area_mean,
         "AspectRatio_mean": aspect_ratio_mean,
         "Donor_mean_intensity_mean": donor_mean,
         "Acceptor_mean_intensity_mean": acceptor_mean,
         "FRET_mean_intensity_mean": fret_mean,
         "FRET_efficiency_mean": fret_eff_mean,
    })

# Create a DataFrame from the list of dictionaries
summary_df = pd.DataFrame(summary_list)

# Save the results to the same folder as "ROI_Data_Summary.csv"
output_file = folder_path / "ROI_Data_Summary.csv"
summary_df.to_csv(output_file, index=False)

print("Summary CSV file has been created at:", output_file)
