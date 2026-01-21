#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute summary statistics from multiple per-file metrics tables and export a single summary CSV.

Assumptions
- Each input file contains the same core columns as the original pipeline output:
  area, Donor_mean_intensity, Acceptor_mean_intensity, FRET_mean_intensity, FRET_efficiency
- In addition, the following columns are present:
  centroid-X, centroid-Y, major_axis_length, minor_axis_length

What this script does
1. Finds all files matching "Summary_metrics_*.csv" in the specified folder and sorts them
   in ascending order by the trailing number in the filename.
2. For each file (using all rows; no filtering), computes:
   - Number of data points (Puncta_Number)
   - Mean of area
   - Mean of inverse aspect ratio (minor_axis_length / major_axis_length)
   - Mean of Donor_mean_intensity
   - Mean of Acceptor_mean_intensity
   - Mean of FRET_mean_intensity
   - Mean of FRET_efficiency
   - Pearson correlation between Acceptor_mean_intensity and FRET_efficiency
     (set to 0 when fewer than two data points are available)
3. Writes the results as a CSV with a header in the specified column order.
"""
import os
import glob
import json
from pathlib import Path
import pandas as pd
import numpy as np

# -----------------------------
# Resolve project root (this script is in project-root/scripts/)
project_root = Path(__file__).resolve().parents[1]

# Load config (sample_id, result_tag)
config_path = project_root / "config" / "fret_proximity_ratio.json"
with config_path.open("r", encoding="utf-8") as f:
    cfg = json.load(f)

sample_id = cfg["sample_id"]
result_tag = cfg["result_tag"]

# Input: project-root/data/csv/{sample_id}
input_folder = project_root / "data" / "csv" / sample_id
if not input_folder.exists():
    raise FileNotFoundError(
        f"Input folder not found: {input_folder}\n"
        f"Check 'sample_id' in config: {config_path}"
    )
if not input_folder.is_dir():
    raise NotADirectoryError(f"Input path is not a directory: {input_folder}")

# Output: project-root/results/summary
output_folder = project_root / "results" / "summary"
output_folder.mkdir(parents=True, exist_ok=True)  # ensure folder exists
# -----------------------------

# Collect the CSV files.
pattern = str(input_folder / "Summary_metrics_*.csv")
csv_files = glob.glob(pattern)
if len(csv_files) == 0:
    raise FileNotFoundError(
        "No input CSV files found.\n"
        f"Expected pattern: Summary_metrics_*.csv\n"
        f"Search folder: {input_folder}"
    )

# Function to sort by the number written at the end of the filename
def extract_number(file_path):
    base = os.path.basename(file_path)
    base_no_ext = base.rsplit(".csv", 1)[0]
    try:
        num = int(base_no_ext.split('_')[-1])
    except ValueError:
        num = float('inf')
    return num

csv_files_sorted = sorted(csv_files, key=extract_number)

# List to store the results
results = []

for csv_path in csv_files_sorted:
    # Extract the identifier (the trailing number) from the filename
    base = os.path.basename(csv_path)
    try:
        file_number = int(base.rsplit(".csv", 1)[0].split('_')[-1])
    except ValueError:
        continue  # Skip if the number cannot be obtained

    # Read the CSV file
    df = pd.read_csv(csv_path)
    df_all = df

    # Number of data points
    puncta_number = len(df_all)
    
    # If no data exist, set everything to 0
    if df_all.empty:
       area_mean = 0
       inv_aspect_mean = 0
       donor_mean = 0
       acceptor_mean = 0
       fret_mean = 0
       fret_eff_mean = 0
       corr = 0
    else:
        # Compute each statistic
        area_mean = df_all['area'].mean()
        inv_aspect = df_all['minor_axis_length'] / df_all['major_axis_length'] # inverse aspect ratio = minor_axis_length / major_axis_length
        inv_aspect_mean = inv_aspect.mean()
        donor_mean = df_all['Donor_mean_intensity'].mean()
        acceptor_mean = df_all['Acceptor_mean_intensity'].mean()
        fret_mean = df_all['FRET_mean_intensity'].mean()
        fret_eff_mean = df_all['FRET_efficiency'].mean()
     
        # Compute the Pearson correlation coefficient (set to 0 if fewer than 2 data points)
        if puncta_number < 2:
            corr = 0
        else:
            corr = np.corrcoef(df_all['Acceptor_mean_intensity'], df_all['FRET_efficiency'])[0, 1]

    # Append the results to the list in the specified order    
    results.append([
        file_number,
        puncta_number,
        area_mean,
        inv_aspect_mean,
        donor_mean,
        acceptor_mean,
        fret_mean,
        fret_eff_mean,
        corr  # Pearson correlation coefficient between Acceptor_mean_intensity and FRET_efficiency
    ])

# Sort the results list in ascending order by the file identifier
results_sorted = sorted(results, key=lambda x: x[0])

# Append the input folder name to the output filename (retrieve the last folder name).
output_csv = output_folder / f"Data_Summary_{result_tag}.csv"

# Output as a CSV file (with headers, 9 columns)
header = [
    "File",
    "Puncta_Number",
    "area_mean",
    "inverse_aspect_mean",
    "Donor_mean_intensity_mean",
    "Acceptor_mean_intensity_mean",
    "FRET_mean_intensity_mean",
    "FRET_efficiency_mean",
    "CC_Acceptor_FRETeff"  # Pearson correlation coefficient
]

df_out = pd.DataFrame(results_sorted, columns=header)
df_out.to_csv(output_csv, index=False)
