import pandas as pd
import matplotlib.pyplot as plt
import os

# ==========================================================
# CONFIGURATION
# ==========================================================

CSV_FILE = "descriptors.csv"
OUTPUT_FOLDER = "plots_rings"
RING_RANGE = list(range(3, 23))  # rings from 3 to 10

# Create output folder automatically
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# ==========================================================
# CSV READING
# ==========================================================

df = pd.read_csv(CSV_FILE)

ring_cols = [f"ring_{n}" for n in RING_RANGE]

# Safety check
for col in ring_cols:
    if col not in df.columns:
        raise ValueError(f"Column {col} not found in the CSV!")

# ==========================================================
# 1) AVERAGE RING DISTRIBUTION (SAVE AS PNG)
# ==========================================================

mean_rings = df[ring_cols].mean().values

plt.figure(figsize=(8, 6))
plt.bar(RING_RANGE, mean_rings)
plt.xlabel("Ring size", fontsize=16, fontweight="bold")
plt.ylabel("Average fraction", fontsize=16, fontweight="bold")
plt.title("Average ring distribution (3–22)", fontsize=18, fontweight="bold")
plt.grid(True)

plt.xticks(fontsize=14, fontweight="bold")
plt.yticks(fontsize=14, fontweight="bold")

for spine in plt.gca().spines.values():
    spine.set_linewidth(1.5)

plt.tight_layout()

mean_path = os.path.join(OUTPUT_FOLDER, "average_ring_distribution.png")
plt.savefig(mean_path, dpi=300)
plt.close()

print("✅ Plot saved:", mean_path)

# ==========================================================
# 2) RING DISTRIBUTION FOR EACH STRUCTURE (SAVE AS PNG)
# ==========================================================

for idx in range(len(df)):
    values = df.loc[idx, ring_cols].values
    name = df.loc[idx, "structure"]

    plt.figure(figsize=(8, 6))
    plt.bar(RING_RANGE, values)
    plt.xlabel("Ring size", fontsize=16, fontweight="bold")
    plt.ylabel("Fraction", fontsize=16, fontweight="bold")
    plt.title(f"Ring distribution – {name}", fontsize=18, fontweight="bold")
    plt.grid(True)

    plt.xticks(fontsize=14, fontweight="bold")
    plt.yticks(fontsize=14, fontweight="bold")

    for spine in plt.gca().spines.values():
        spine.set_linewidth(1.5)

    plt.tight_layout()

    clean_name = os.path.splitext(name)[0]
    out_path = os.path.join(
        OUTPUT_FOLDER,
        f"rings_{idx:04d}_{clean_name}.png"
    )

    plt.savefig(out_path, dpi=200)
    plt.close()

print("✅ Individual structure plots saved in folder:", OUTPUT_FOLDER)

# ==========================================================
# 3) HEATMAP (SAVE AS PNG)
# ==========================================================

plt.figure(figsize=(10, 6))
plt.imshow(df[ring_cols].values, aspect="auto")
plt.xlabel("Ring type (3–10)", fontsize=16, fontweight="bold")
plt.ylabel("Structure", fontsize=16, fontweight="bold")
plt.xticks(range(len(RING_RANGE)), RING_RANGE, fontsize=14, fontweight="bold")
plt.yticks(fontsize=14, fontweight="bold")
plt.title("Heatmap of ring distribution", fontsize=18, fontweight="bold")

cbar = plt.colorbar()
cbar.set_label("Fraction", fontsize=14, fontweight="bold")
cbar.ax.tick_params(labelsize=12)
for tick in cbar.ax.get_yticklabels():
    tick.set_fontweight("bold")

for spine in plt.gca().spines.values():
    spine.set_linewidth(1.5)

plt.tight_layout()

heatmap_path = os.path.join(OUTPUT_FOLDER, "ring_distribution_heatmap.png")
plt.savefig(heatmap_path, dpi=300)
plt.close()

print("✅ Heatmap saved:", heatmap_path)

# ==========================================================
# 4) QUICK TEXT REPORT (OPTIONAL)
# ==========================================================

report_path = os.path.join(OUTPUT_FOLDER, "ring_distribution_summary.txt")

with open(report_path, "w") as f:
    f.write("Statistical summary of ring distribution (3–10)\n\n")
    f.write(df[ring_cols].describe().to_string())

print("✅ Statistical report saved:", report_path)

print("\n🎉 ALL PLOTS WERE GENERATED SUCCESSFULLY!")
