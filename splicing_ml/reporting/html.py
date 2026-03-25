from __future__ import annotations

"""HTML report renderer for subset configurations.

The main entry point is ``write_subset_html_report``, which writes a
self-contained HTML file (plus a sibling JSON payload file) for one subset
configuration.

Key improvement over the original monolithic f-string: the eight repeated
``Plotly.newPlot(...)`` call patterns are extracted into named JavaScript
helper functions (``_PLOTLY_JS_HELPERS``), assembled as a plain string rather
than an f-string so that curly braces inside the JS do not need escaping.
This reduces the renderReports body by ~80 lines and makes each plot type
independently auditable.
"""

import html
import json
from pathlib import Path
from typing import Any

from ..utils import safe_json, vlog
from .payload import build_task_plot_payload

__all__ = [
    "config_key",
    "slugify_config_key",
    "write_subset_html_report",
    "generate_html_reports",
]

# ---------------------------------------------------------------------------
# Config key helpers
# ---------------------------------------------------------------------------


def config_key(result: dict[str, Any]) -> tuple[str, str, str, str]:
    """Extract the (event_type, transcript_filter, variability, group_col) key."""
    cfg = result.get("config", {})
    return (
        str(cfg.get("event_type", "NA")),
        str(cfg.get("transcript_filter", "NA")),
        str(cfg.get("variability", "NA")),
        str(cfg.get("group_col", "NA")),
    )


def slugify_config_key(key: tuple[str, str, str, str]) -> str:
    """Return a filename-safe slug from a subset config key."""
    raw = "_".join(key)
    return "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in raw)


# ---------------------------------------------------------------------------
# JavaScript helper functions (plain string — {} in JS needs no escaping here)
# ---------------------------------------------------------------------------

_PLOTLY_JS_HELPERS = """
    // ── Metric heatmap ──────────────────────────────────────────────────────
    function plotMetricHeatmap(divId, hm, title) {
        const el = document.getElementById(divId);
        if (!el) return;
        if ((hm.models || []).length === 0 || (hm.metrics || []).length === 0) {
            el.innerHTML = '<p>No heatmap data available.</p>';
            return;
        }
        Plotly.newPlot(
            divId,
            [{
                type: 'heatmap',
                x: hm.metrics,
                y: hm.models,
                z: hm.z,
                customdata: hm.sd,
                text: hm.text,
                colorscale: 'Viridis',
                texttemplate: '%{text}',
                hovertemplate: 'Model=%{y}<br>Metric=%{x}<br>Mean=%{z:.4f}<br>SD=%{customdata:.4f}<extra></extra>'
            }],
            {
                title: title,
                xaxis: { title: 'Metric' },
                yaxis: { title: 'Model', automargin: true }
            },
            { responsive: true }
        );
    }

    // ── Observed vs Predicted scatter ────────────────────────────────────────
    function plotScatterOvP(divId, xData, yData, modelName, xLabel, yLabel, title) {
        if (!xData.length || !yData.length) return;
        const minVal = Math.min(...xData, ...yData);
        const maxVal = Math.max(...xData, ...yData);
        Plotly.newPlot(
            divId,
            [
                {
                    x: xData, y: yData, mode: 'markers', type: 'scatter',
                    name: modelName, marker: { opacity: 0.5, size: 5 }
                },
                {
                    x: [minVal, maxVal], y: [minVal, maxVal], mode: 'lines',
                    type: 'scatter', name: 'y=x',
                    line: { color: '#666', width: 1, dash: 'dot' }, hoverinfo: 'skip'
                }
            ],
            {
                title: { text: title, font: { size: 13 } },
                margin: { l: 60, r: 20, t: 45, b: 55 },
                xaxis: { title: xLabel, automargin: true },
                yaxis: { title: yLabel, scaleanchor: 'x', scaleratio: 1, automargin: true }
            },
            { responsive: true }
        );
    }

    // ── PSI distribution histogram ───────────────────────────────────────────
    function plotPsiHistogram(divId, psiVals, title, shapes) {
        Plotly.newPlot(
            divId,
            [{ x: psiVals, type: 'histogram', nbinsx: 60 }],
            {
                title: title,
                xaxis: { title: 'PSI' },
                yaxis: { title: 'Count' },
                shapes: shapes || []
            },
            { responsive: true }
        );
    }

    // ── Confusion matrix heatmap ─────────────────────────────────────────────
    function plotConfusionMatrix(divId, cm, modelName) {
        Plotly.newPlot(
            divId,
            [{
                z: cm,
                x: ['Pred 0', 'Pred 1'],
                y: ['True 0', 'True 1'],
                type: 'heatmap',
                colorscale: 'Blues',
                showscale: true,
                text: cm,
                texttemplate: '%{text}'
            }],
            { title: 'Confusion matrix - model=' + modelName },
            { responsive: true }
        );
    }
"""

# ---------------------------------------------------------------------------
# HTML report writer
# ---------------------------------------------------------------------------


def write_subset_html_report(
    config: tuple[str, str, str, str],
    task_results: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write one standalone HTML report for a subset configuration.

    A sibling ``<name>.json`` file containing the serialised plotting payload
    is written alongside the HTML; the JS fetches it at runtime so the HTML
    itself remains small.
    """
    event_type, transcript_filter, variability, group_col = config
    title = (
        f"Subset Report: event_type={event_type}, "
        f"transcript_filter={transcript_filter}, variability={variability}, "
        f"group_col={group_col}"
    )
    payload = [build_task_plot_payload(tr) for tr in task_results]
    payload_path = output_path.with_suffix(".json")
    payload_path.write_text(json.dumps(safe_json(payload)), encoding="utf-8")
    payload_file = payload_path.name

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(title)}</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 20px; }}
    .card {{ border: 1px solid #ddd; border-radius: 8px; padding: 14px; margin: 14px 0; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #ddd; padding: 6px; text-align: left; font-size: 0.9rem; }}
    th {{ background: #f6f6f6; }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <div id="reports"></div>
  <script>
{_PLOTLY_JS_HELPERS}

    const payloadFile = {json.dumps(payload_file)};
    const container = document.getElementById('reports');

    function createParamsTable(rows) {{
      if (!rows || rows.length === 0) return '<p>No tuning rows available.</p>';
      let out = '<table><thead><tr><th>Model</th><th>Selected in folds</th><th>Most frequent best params</th><th>Frequency</th></tr></thead><tbody>';
      for (const r of rows) {{
        out += `<tr><td>${{r.model_name}}</td><td>${{r.selected_in_folds}}</td><td><code>${{r.most_frequent_best_params}}</code></td><td>${{r.frequency}}</td></tr>`;
      }}
      return out + '</tbody></table>';
    }}

    async function loadPayload(path) {{
      const resp = await fetch(path, {{ cache: 'no-store' }});
      if (!resp.ok) throw new Error(`Failed to load JSON payload: ${{resp.status}}`);
      return await resp.json();
    }}

    function renderReports(reportData) {{
      for (const tr of reportData) {{
        const card = document.createElement('div');
        card.className = 'card';
        card.innerHTML = `<h2>Task: ${{tr.task}} (status=${{tr.status}})</h2><p>Primary metric: ${{tr.primary_metric}} | n_samples: ${{tr.n_samples}}</p>`;

        const p0 = document.createElement('div'); p0.id = `p0_${{tr.task}}`; p0.style.height = '420px';
        const p1 = document.createElement('div'); p1.id = `p1_${{tr.task}}`; p1.style.height = '320px';
        const p2 = document.createElement('div'); p2.id = `p2_${{tr.task}}`; p2.style.height = 'auto';
        const p4 = document.createElement('div'); p4.id = `p4_${{tr.task}}`; p4.style.height = '360px';
        card.appendChild(p0); card.appendChild(p1); card.appendChild(p2); card.appendChild(p4);

        const h3 = document.createElement('h3'); h3.textContent = 'Important parameters'; card.appendChild(h3);
        const ptab = document.createElement('div'); ptab.innerHTML = createParamsTable(tr.important_params); card.appendChild(ptab);
        container.appendChild(card);

        // ── Metric heatmap(s) ────────────────────────────────────────────────
        if (tr.task === 'regression') {{
          p0.style.height = '760px';
          const hmOrig  = tr.metric_heatmap_original || {{ models: [], metrics: [], z: [], sd: [], text: [] }};
          const hmLogit = tr.metric_heatmap_logit    || {{ models: [], metrics: [], z: [], sd: [], text: [] }};
          p0.innerHTML = '';
          const hmOrigDiv  = document.createElement('div'); hmOrigDiv.id  = `p0_orig_${{tr.task}}`;  hmOrigDiv.style.height  = '360px';
          const hmLogitDiv = document.createElement('div'); hmLogitDiv.id = `p0_logit_${{tr.task}}`; hmLogitDiv.style.height = '360px';
          p0.appendChild(hmOrigDiv); p0.appendChild(hmLogitDiv);
          plotMetricHeatmap(hmOrigDiv.id,  hmOrig,  'Regression metrics heatmap (original PSI scale)');
          plotMetricHeatmap(hmLogitDiv.id, hmLogit, 'Regression metrics heatmap (logit PSI scale)');
        }} else {{
          p0.style.height = '420px';
          const heatmap = tr.metric_heatmap || {{ models: [], metrics: [], z: [], sd: [], text: [] }};
          plotMetricHeatmap(p0.id, heatmap, 'All metrics heatmap (mean +/- sd over folds)');
        }}

        // ── Model-wise fold distribution ─────────────────────────────────────
        const boxTraces = [];
        for (const modelName of tr.model_order) {{
          const yVals = tr.model_fold_points
            .filter(p => p.model_name === modelName)
            .map(p => p.primary_score);
          boxTraces.push({{
            type: 'box', name: modelName, y: yVals,
            boxpoints: 'all', jitter: 0.35, pointpos: 0,
            marker: {{ size: 6, opacity: 0.75 }}, line: {{ width: 1 }}
          }});
        }}
        const foldLineTraces = tr.fold_lines
          .filter(fl => fl.x.length >= 2)
          .map(fl => ({{
            x: fl.x, y: fl.y, mode: 'lines+markers', type: 'scatter',
            name: `Fold ${{fl.outer_fold}}`, showlegend: false,
            marker: {{ size: 5, opacity: 0.65 }},
            line: {{ width: 1, color: '#666' }},
            hovertemplate: 'Fold ' + fl.outer_fold + '<br>%{{x}}: %{{y}}<extra></extra>'
          }}));
        Plotly.newPlot(
          p1.id, [...boxTraces, ...foldLineTraces],
          {{
            title: `Model-wise distribution of ${{tr.primary_metric}}`,
            xaxis: {{ title: 'Model type' }},
            yaxis: {{ title: tr.primary_metric }}
          }},
          {{ responsive: true }}
        );

        // ── Predictions / distributions ──────────────────────────────────────
        if (tr.task === 'regression') {{
          const byModelScale = tr.prediction_by_model_scale || {{}};
          p2.innerHTML = '';
          const regHeader = document.createElement('h3');
          regHeader.textContent = 'Observed vs Predicted faceted by model and scale';
          p2.appendChild(regHeader);
          const regGrid = document.createElement('div');
          regGrid.style.display = 'grid';
          regGrid.style.gridTemplateColumns = `repeat(${{Math.max(1, tr.model_order.length)}}, minmax(260px, 1fr))`;
          regGrid.style.gap = '12px';
          regGrid.style.alignItems = 'flex-start';
          p2.appendChild(regGrid);

          const scales = [
            {{ key: 'original', xLabel: 'Observed PSI',         yLabel: 'Predicted PSI' }},
            {{ key: 'logit',    xLabel: 'Observed logit(PSI)',   yLabel: 'Predicted logit(PSI)' }}
          ];
          let hasRegressionPanels = false;
          for (const s of scales) {{
            for (const modelName of tr.model_order) {{
              const m  = byModelScale[modelName] || {{ original: {{ y_true: [], y_pred: [] }}, logit: {{ y_true: [], y_pred: [] }} }};
              const sx = (m[s.key] || {{}}).y_true || [];
              const sy = (m[s.key] || {{}}).y_pred || [];
              if (!sx.length || !sy.length) continue;
              hasRegressionPanels = true;
              const regDiv = document.createElement('div');
              const safeModelName = String(modelName).replace(/[^a-zA-Z0-9_]/g, '_');
              regDiv.id = `ovp_${{tr.task}}_${{safeModelName}}_${{s.key}}`;
              regDiv.style.minWidth = '0'; regDiv.style.height = '360px';
              regGrid.appendChild(regDiv);
              plotScatterOvP(
                regDiv.id, sx, sy, modelName, s.xLabel, s.yLabel,
                `Observed vs Predicted - ${{modelName}} (${{s.key}})`
              );
            }}
          }}
          if (!hasRegressionPanels) {{
            const noData = document.createElement('p');
            noData.textContent = 'No observed/predicted data available for regression faceting.';
            p2.appendChild(noData);
          }}

          const dist    = tr.response_distribution || {{}};
          const psiVals = dist.psi_sample || [];
          const thr = dist.binarization_thresholds || [];
          const regShapes = (thr.length === 2) ? [
            {{ type: 'rect', x0: 0, x1: thr[0], y0: 0, y1: 1, yref: 'paper',
               fillcolor: 'rgba(180,0,0,0.10)', line: {{ width: 0 }},
               layer: 'below' }},
            {{ type: 'rect', x0: thr[1], x1: 1, y0: 0, y1: 1, yref: 'paper',
               fillcolor: 'rgba(180,0,0,0.10)', line: {{ width: 0 }},
               layer: 'below' }}
          ] : [];
          plotPsiHistogram(p4.id, psiVals, 'Response (PSI) distribution (shaded = excluded extremes)', regShapes);

        }} else {{
          // Classification: predicted probability histograms + threshold lines.
          const byModel      = tr.prediction_by_model || {{}};
          const thresholdRows = tr.threshold_by_model || [];
          const palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'];
          const colorByModel = {{}};
          const probTraces   = [];

          for (const modelName of tr.model_order) {{
            colorByModel[modelName] = palette[Object.keys(colorByModel).length % palette.length];
            const m = byModel[modelName] || {{ y_pred: [] }};
            if (!(m.y_pred || []).length) continue;
            probTraces.push({{
              x: m.y_pred, type: 'histogram', nbinsx: 40,
              name: modelName, opacity: 0.5,
              marker: {{ color: colorByModel[modelName] }}
            }});
          }}

          const thresholdShapes = [];
          const thresholdAnn    = [];
          for (let i = 0; i < thresholdRows.length; i++) {{
            const row = thresholdRows[i];
            const thr = Number(row.mean_threshold);
            if (Number.isFinite(thr)) {{
              const modelName = String(row.model_name || 'model');
              const c = colorByModel[modelName] || '#111';
              thresholdShapes.push({{ type: 'line', x0: thr, x1: thr, y0: 0, y1: 1, yref: 'paper', line: {{ color: c, width: 2, dash: 'dot' }} }});
              thresholdAnn.push({{
                x: thr, y: 1.02 + (0.045 * (i % 2)), yref: 'paper',
                text: `${{modelName}} thr=${{thr.toFixed(3)}}`,
                showarrow: false, xanchor: 'left', font: {{ size: 10, color: c }}
              }});
            }}
          }}
          Plotly.newPlot(
            p2.id, probTraces,
            {{
              title: 'Predicted probability distribution by model (model thresholds shown)',
              xaxis: {{ title: 'P(class=1)' }}, yaxis: {{ title: 'Count' }},
              barmode: 'overlay', shapes: thresholdShapes, annotations: thresholdAnn
            }},
            {{ responsive: true }}
          );

          // Per-model confusion matrices.
          const cmHeader = document.createElement('h3');
          cmHeader.textContent = 'Confusion matrices by model'; card.appendChild(cmHeader);
          const confWrap = document.createElement('div'); card.appendChild(confWrap);
          const confByModel = tr.confusion_by_model || {{}};
          for (const modelName of Object.keys(confByModel)) {{
            const cmDiv = document.createElement('div');
            cmDiv.id = `cm_${{tr.task}}_${{modelName}}`; cmDiv.style.height = '300px';
            confWrap.appendChild(cmDiv);
            plotConfusionMatrix(cmDiv.id, confByModel[modelName], modelName);
          }}

          // PSI distribution with binarization threshold markers.
          const dist = tr.response_distribution || {{}};
          const psiVals = dist.psi_sample || [];
          const thr  = dist.binarization_thresholds || [];
          const shapes = (thr.length === 2) ? [
            {{ type: 'rect', x0: thr[0], x1: thr[1], y0: 0, y1: 1, yref: 'paper',
               fillcolor: 'rgba(180,0,0,0.10)', line: {{ width: 0 }},
               layer: 'below' }},
            {{ type: 'line', x0: thr[0], x1: thr[0], y0: 0, y1: 1, yref: 'paper', line: {{ color: '#c00', width: 2, dash: 'dash' }} }},
            {{ type: 'line', x0: thr[1], x1: thr[1], y0: 0, y1: 1, yref: 'paper', line: {{ color: '#c00', width: 2, dash: 'dash' }} }}
          ] : [];
          plotPsiHistogram(p4.id, psiVals, 'Response (PSI) distribution (shaded = excluded mid-range)', shapes);
        }}
      }}
    }}

    loadPayload(payloadFile)
      .then(renderReports)
      .catch((err) => {{
        const msg = document.createElement('div');
        msg.className = 'card';
        msg.innerHTML = `<h2>Report Load Error</h2><p>${{String(err)}}</p><p>If you opened this HTML via file://, serve the folder over HTTP so fetch() can read the JSON payload.</p>`;
        container.appendChild(msg);
      }});
  </script>
</body>
</html>
"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html_text, encoding="utf-8")


# ---------------------------------------------------------------------------
# Batch report generator
# ---------------------------------------------------------------------------


def generate_html_reports(
    all_results: list[dict[str, Any]], out_dir: Path, verbose: bool = False
) -> None:
    """Generate one HTML report per subset configuration."""
    import collections

    grouped: dict[tuple[str, str, str, str], list[dict[str, Any]]] = (
        collections.defaultdict(list)
    )
    for result in all_results:
        grouped[config_key(result)].append(result)

    reports_dir = out_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    for key, task_results in grouped.items():
        report_path = reports_dir / f"subset_report_{slugify_config_key(key)}.html"
        write_subset_html_report(key, task_results, report_path)
        vlog(verbose, f"Wrote HTML report: {report_path}", level="info")
