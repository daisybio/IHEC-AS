"""Push the Fig2B-style comparison data to W&B as a native Table + bar charts.

Called from 07-2-ml-global-comparison.Rmd via system2(). Logs the raw tidy data and
lets wandb build the charts itself (interactive/filterable in the UI) rather than
uploading a rendered image -- avoids the PIL dependency wandb.Image() needs (not
installed in this project's env) and needs no pandas/polars: stdlib csv is enough
for a table this small.
"""
import csv
import sys
from collections import defaultdict

import wandb

csv_path = sys.argv[1]

with open(csv_path, newline="") as f:
    rows = list(csv.DictReader(f))

columns = list(rows[0].keys())
data = [[r[c] for c in columns] for r in rows]

run = wandb.init(project="splicing-ml", name="fig2b-style-comparison", job_type="comparison-figure", reinit=True)

table = wandb.Table(columns=columns, data=data)
run.log({"comparison_table": table})

# One bar chart per (Event Type, Variability, Metric). wandb.plot.bar has no native
# grouped-bar support (one bar per table row), so the x-axis label combines the
# held-out condition and feature set, e.g. "Chr | Epigenetic" vs "Cell | All".
groups = defaultdict(list)
for r in rows:
    groups[(r["Event Type"], r["Variability"], r["Metric"])].append(r)

for (event_type, variability, metric), group_rows in groups.items():
    chart_data = [[f"{r['Test on']} | {r['Features']}", float(r["Value"])] for r in group_rows]
    chart_table = wandb.Table(columns=["label", "Value"], data=chart_data)
    panel_key = f"{event_type}_{variability}/{metric}"
    run.log({panel_key: wandb.plot.bar(chart_table, "label", "Value", title=f"{event_type} {variability} {metric}")})

run.finish()
