import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon
import starbars

def significance_stars(p_value):
    if p_value <= 0.00001:
        return "*****"
    elif p_value <= 0.0001:
        return "****"
    elif p_value <= 0.001:
        return "***"
    elif p_value <= 0.01:
        return "**"
    elif p_value <= 0.05:
        return "*"
    else:
        return f"{p_value:.3f}"  # Not significant

def draw_boxplot(combined_df, save_path):
    labels = ["Raw", "eHOMD", "HROM"]

    plt.figure(figsize=(8, 8))

    custom_palette = {
        labels[0]: "white",
        labels[1]: "#1f77b4",
        labels[2]: "#ff7f0e"
    }

    # Plot boxplots
    ax = sns.boxplot(
        data=combined_df,
        x="Metric",
        y="Value",
        hue="Group",
        palette=custom_palette
    )

    plt.title("")
    plt.xlabel("Metrics")
    plt.ylabel("Value")
    plt.xticks(rotation=45, ha="right")

    # Add statistical annotations using Wilcoxon signed-rank test
    annotations = []
    for metric in ["F1 Score", "Precision", "Recall"]:
        for group_pair in [(labels[0], labels[1]), (labels[0], labels[2]), (labels[1], labels[2])]:  
            group1 = combined_df[(combined_df["Metric"] == metric) & (combined_df["Group"] == group_pair[0])]["Value"]
            group2 = combined_df[(combined_df["Metric"] == metric) & (combined_df["Group"] == group_pair[1])]["Value"]

            # Perform Wilcoxon signed-rank test
            if len(group1) > 0 and len(group2) > 0 and len(group1) == len(group2):
                stat, p_value = wilcoxon(group1, group2, alternative="less")
                annotations.append(((group_pair[0], metric), (group_pair[1], metric), significance_stars(p_value)))

    # Adding statistical annotations
    starbars.draw_annotation(annotations)

    plt.tight_layout()
    
    # Save the plot instead of showing it
    plt.savefig(save_path, dpi=300)
    plt.close()  # Close the figure to free memory

# Define file paths
file_paths = {
    "common_SNP": "/home/zunuan/for_JHC/validation/common_SNP_dataframe.tsv",
    "common_INDEL": "/home/zunuan/for_JHC/validation/common_INDEL_dataframe.tsv",
    "rare_SNP": "/home/zunuan/for_JHC/validation/rare_SNP_dataframe.tsv",
    "rare_INDEL": "/home/zunuan/for_JHC/validation/rare_INDEL_dataframe.tsv"
}

# Loop through files and generate boxplots
for key, file_path in file_paths.items():
    combined_df = pd.read_csv(file_path, sep="\t")
    save_path = f"{key}_boxplot.png"  # Save as PNG
    draw_boxplot(combined_df, save_path)

    print(f"Saved boxplot to {save_path}")
