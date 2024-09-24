""" Observations
Capomulin and Ramicane show the lowest average tumor volumes among the treatment regimens, 
suggesting these drugs may be the most effective in reducing tumor growth.

Infubinol shows an outlier in the final tumor volume distribution, indicating that while it may be effective in some cases, 
there are instances where it fails to control tumor growth effectively.

A strong positive correlation (0.84) is observed between mouse weight and tumor volume for the Capomulin regimen, 
suggesting that heavier mice tend to have larger tumor volumes during the treatment.
"""
# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Load datasets
mouse_metadata = pd.read_csv('/Users/ryanpope/Pymaceuticals/data/Mouse_metadata.csv')
study_results = pd.read_csv('/Users/ryanpope/Pymaceuticals/data/Study_results.csv')

# Merge datasets on "Mouse ID"
merged_data = pd.merge(study_results, mouse_metadata, on="Mouse ID")

# Remove duplicates and clean data
cleaned_data = merged_data.drop_duplicates(subset=["Mouse ID", "Timepoint"])

# Generate summary statistics
summary_stats = cleaned_data.groupby("Drug Regimen").agg(
    mean_tumor_volume=("Tumor Volume (mm3)", "mean"),
    median_tumor_volume=("Tumor Volume (mm3)", "median"),
    tumor_volume_variance=("Tumor Volume (mm3)", "var"),
    tumor_volume_std=("Tumor Volume (mm3)", "std"),
    tumor_volume_sem=("Tumor Volume (mm3)", "sem")
)

# Bar charts for total timepoints by drug regimen
timepoint_counts = cleaned_data["Drug Regimen"].value_counts()
timepoint_counts.plot(kind='bar', figsize=(10,6), title='Total Timepoints for Each Drug Regimen (Pandas)')
plt.ylabel('Number of Timepoints')
plt.show()

plt.figure(figsize=(10,6))
plt.bar(timepoint_counts.index, timepoint_counts.values)
plt.title('Total Timepoints for Each Drug Regimen (Matplotlib)')
plt.xlabel('Drug Regimen')
plt.ylabel('Number of Timepoints')
plt.xticks(rotation=45)
plt.show()

# Pie charts for sex distribution
sex_counts = cleaned_data["Sex"].value_counts()
sex_counts.plot(kind='pie', autopct='%1.1f%%', figsize=(6,6), title='Distribution of Male vs Female Mice (Pandas)')
plt.show()

plt.figure(figsize=(6,6))
plt.pie(sex_counts, labels=sex_counts.index, autopct='%1.1f%%')
plt.title('Distribution of Male vs Female Mice (Matplotlib)')
plt.show()

# Quartiles, IQR, and box plots for final tumor volumes of promising treatments
promising_regimens = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]
max_timepoints = cleaned_data.groupby("Mouse ID")["Timepoint"].max().reset_index()
final_tumor_volume_df = pd.merge(max_timepoints, cleaned_data, on=["Mouse ID", "Timepoint"])

final_tumor_volumes = final_tumor_volume_df[final_tumor_volume_df["Drug Regimen"].isin(promising_regimens)]

plt.figure(figsize=(10,6))
final_tumor_volumes.boxplot(column="Tumor Volume (mm3)", by="Drug Regimen", grid=False, patch_artist=True)
plt.title('Final Tumor Volume Distribution by Treatment')
plt.suptitle('')
plt.ylabel('Tumor Volume (mm3)')
plt.show()

# Line plot for a single mouse treated with Capomulin
capomulin_mouse = cleaned_data[(cleaned_data["Drug Regimen"] == "Capomulin") & (cleaned_data["Mouse ID"] == "l509")]
plt.figure(figsize=(8,6))
plt.plot(capomulin_mouse["Timepoint"], capomulin_mouse["Tumor Volume (mm3)"], marker="o")
plt.title("Tumor Volume vs. Timepoint for Mouse l509 (Capomulin)")
plt.xlabel("Timepoint (days)")
plt.ylabel("Tumor Volume (mm3)")
plt.show()

# Scatter plot and linear regression for Capomulin regimen
capomulin_grouped = cleaned_data[cleaned_data["Drug Regimen"] == "Capomulin"].groupby("Mouse ID").mean()

weight = capomulin_grouped["Weight (g)"]
tumor_volume = capomulin_grouped["Tumor Volume (mm3)"]
slope, intercept, r_value, p_value, std_err = linregress(weight, tumor_volume)

regression_line = slope * weight + intercept

plt.figure(figsize=(8,6))
plt.scatter(weight, tumor_volume)
plt.plot(weight, regression_line, color="red")
plt.title("Mouse Weight vs. Average Tumor Volume (Capomulin) with Regression Line")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()

# Output correlation coefficient
print(f"Correlation Coefficient: {r_value}")
