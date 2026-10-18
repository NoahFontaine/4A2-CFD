from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


csv_path = Path(__file__).with_name("edge_vector_errors.csv")
data = pd.read_csv(csv_path)

numeric_data = data.select_dtypes(include="number")
if numeric_data.empty:
	raise ValueError(f"No numeric columns found in {csv_path.name}")

fig, axes = plt.subplots(
	len(numeric_data.columns), 1,
	figsize=(10, max(6, 2.5 * len(numeric_data.columns))),
	sharex=True,
	squeeze=False,
)
for ax, column in zip(axes[:, 0], numeric_data.columns):
	numeric_data[column].plot(ax=ax, marker=".", legend=False)
	ax.set_ylabel(str(column))
	ax.grid(True, alpha=0.3)

axes[-1, 0].set_xlabel("Sample")
fig.suptitle("Edge Vector Errors")
plt.tight_layout()
plt.savefig("edge_vector_errors.png")
