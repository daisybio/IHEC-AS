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
from datetime import datetime, timezone
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
# Model card
# ---------------------------------------------------------------------------

# Caveat predicates: each takes (group_col, model_names) and, if it applies,
# contributes one <li> to the model card. Kept in sync with CLAUDE.md's
# "Known Pitfalls" section by hand -- there is no single source of truth to
# derive this from automatically, so if a pitfall there changes, update here.
_MODEL_CARD_CAVEATS: list[Any] = [
    (
        lambda group_col, models: group_col == "ontology",
        'group_col="ontology" outer CV leaks event identity for high-capacity '
        "models (xgb/lgbm/tabicl AUROC 0.97–1.0 on both Low and High "
        "variability) — rows from the same event scatter across ontology "
        "labels, so the model can fingerprint the event instead of "
        "generalizing to a genuinely unseen cell type. Hierarchical ontology "
        "supergrouping is applied but confirmed insufficient. Do not cite "
        "tree/tabicl AUROC from this report as real generalization "
        "performance; linear's numbers are unaffected.",
    ),
    (
        lambda group_col, models: {"xgb", "lgbm"} <= models,
        "lgbm vs xgb SHAP: the two agree strongly on feature <em>ranking</em> "
        "(Spearman ρ≈0.885) but lgbm's absolute SHAP magnitude runs "
        "~1.6–35× larger — compare rankings across these two "
        "model families, never raw SHAP magnitudes.",
    ),
    (
        lambda group_col, models: "tabicl" in models,
        "tabicl has no SHAP/feature-attribution support (pretrained, no inner "
        "CV, never routed through the tree-model SHAP path) — expect no "
        "per-feature explanation for it in this report.",
    ),
    (
        lambda group_col, models: "rf" in models,
        '"rf" is XGBoost’s XGBRFClassifier/XGBRFRegressor, not sklearn’s '
        "RandomForest — shares xgb's GPU setup but gets none of the "
        "early-stopping machinery.",
    ),
]


def _build_model_card_html(
    config: tuple[str, str, str, str],
    payload: list[dict[str, Any]],
    run_datetime_display: str,
) -> str:
    """Return the "Model Card" HTML block shown at the top of the report.

    Summarises the subset config and flags known caveats (CV leakage,
    cross-model SHAP comparability, etc.) that apply to the specific
    combination of group_col/models present in this report -- so a reader
    doesn't need CLAUDE.md open to know a number shouldn't be trusted at
    face value.
    """
    event_type, transcript_filter, variability, group_col = config
    model_names: set[str] = set()
    for p in payload:
        model_names.update(p.get("model_order") or [])
    models_display = ", ".join(sorted(model_names)) if model_names else "N/A"

    caveats = [msg for predicate, msg in _MODEL_CARD_CAVEATS if predicate(group_col, model_names)]
    if caveats:
        caveats_html = "<ul style='margin:6px 0 0;padding-left:20px;'>" + "".join(
            f"<li>{c}</li>" for c in caveats
        ) + "</ul>"
    else:
        caveats_html = "<p style='color:#666;margin:6px 0 0;'>No known caveats apply to this subset/model combination.</p>"

    return f"""
  <div class="card" style="background:#fffbea;border-color:#e8d9a0;">
    <h2 style="margin-top:0;">Model Card</h2>
    <table style="width:auto;">
      <tr><th>Event type</th><td>{html.escape(event_type)}</td></tr>
      <tr><th>Transcript filter</th><td>{html.escape(transcript_filter)}</td></tr>
      <tr><th>Variability</th><td>{html.escape(variability)}</td></tr>
      <tr><th>Group column (outer CV)</th><td>{html.escape(group_col)}</td></tr>
      <tr><th>Models</th><td>{html.escape(models_display)}</td></tr>
      <tr><th>Results generated</th><td>{run_datetime_display}</td></tr>
    </table>
    <h3 style="margin-bottom:4px;">Known caveats for this configuration</h3>
    {caveats_html}
  </div>
"""


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
        const nModels = hm.models.length;
        const nMetrics = hm.metrics.length;
        // hm.z/sd/text arrive as model-rows x metric-cols; transpose so metrics
        // become rows (y-axis) and models become columns (x-axis).
        const zT = [], sdT = [], textT = [];
        for (let mi = 0; mi < nMetrics; mi++) {
            const zRow = [], sdRow = [], textRow = [];
            for (let ri = 0; ri < nModels; ri++) {
                zRow.push(hm.z[ri][mi]);
                sdRow.push(hm.sd[ri][mi]);
                textRow.push(hm.text[ri][mi]);
            }
            zT.push(zRow); sdT.push(sdRow); textT.push(textRow);
        }
        // Metrics live on unrelated scales (e.g. RMSE vs AUROC), so color each
        // row (metric) independently via min-max scaling; the actual mean/sd
        // stay available as displayed text and in hover via customdata.
        // Some metrics are lower-is-better (RMSE, MAD) and some are
        // higher-is-better (R2, CCC, AUROC, ...) -- without accounting for
        // that, "best" would map to opposite color ends on different rows.
        // hm.directions[i] (true = lower-is-better) flips the per-row scale
        // so yellow always means "best" and purple "worst". Rows also arrive
        // pre-grouped by direction (server-side sort), so a single divider
        // line cleanly separates the two groups.
        const directions = hm.directions || hm.metrics.map(() => false);
        const yLabels = hm.metrics.map((m, i) => m + (directions[i] ? ' ↓ lower better' : ' ↑ higher better'));
        const zColor = zT.map((row, ri) => {
            const finite = row.filter((v) => v !== null && v !== undefined && !Number.isNaN(v));
            if (!finite.length) return row.map(() => null);
            const lo = Math.min(...finite), hi = Math.max(...finite);
            const lowerIsBetter = !!directions[ri];
            if (hi === lo) return row.map((v) => (v === null || v === undefined ? null : 0.5));
            return row.map((v) => {
                if (v === null || v === undefined) return null;
                const frac = (v - lo) / (hi - lo);
                return lowerIsBetter ? 1 - frac : frac;
            });
        });
        const customdata = zT.map((row, ri) => row.map((v, ci) => [v, sdT[ri][ci]]));

        const nLowerIsBetter = directions.filter(Boolean).length;
        const groupShapes = [];
        if (nLowerIsBetter > 0 && nLowerIsBetter < nMetrics) {
            // Lower-is-better rows are grouped first (bottom of the
            // categorical y-axis), so the boundary sits right above them.
            groupShapes.push({
                type: 'line', xref: 'paper', x0: 0, x1: 1,
                yref: 'y', y0: nLowerIsBetter - 0.5, y1: nLowerIsBetter - 0.5,
                line: { color: '#888', width: 1, dash: 'dot' }
            });
        }

        Plotly.newPlot(
            divId,
            [{
                type: 'heatmap',
                x: hm.models,
                y: yLabels,
                z: zColor,
                zmin: 0,
                zmax: 1,
                customdata: customdata,
                text: textT,
                colorscale: 'Viridis',
                texttemplate: '%{text}',
                showscale: true,
                colorbar: { title: 'Row-scaled:<br>1=best, 0=worst' },
                hovertemplate: 'Model=%{x}<br>Metric=%{y}<br>Mean=%{customdata[0]:.4f}<br>SD=%{customdata[1]:.4f}<extra></extra>'
            }],
            {
                title: title,
                xaxis: { title: 'Model', automargin: true },
                yaxis: { title: 'Metric', automargin: true },
                shapes: groupShapes
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

    // ── Transformed (post-preprocessing) feature distribution ───────────────
    function plotFeatureHistogram(divId, values, featureName) {
        Plotly.newPlot(
            divId,
            [{ x: values, type: 'histogram', nbinsx: 40 }],
            {
                title: { text: featureName, font: { size: 12 } },
                margin: { l: 45, r: 15, t: 35, b: 35 },
                xaxis: { title: 'Transformed value' },
                yaxis: { title: 'Count' }
            },
            { responsive: true }
        );
    }

    // ── ROC + PR combined panel (shared legend) ────────────────────────────
    function plotRocPrCurves(divId, rocByModel, prByModel, prevalence, palette) {
        const el = document.getElementById(divId);
        if (!el) return;
      const roc = rocByModel || {};
      const pr = prByModel || {};
      const BASELINE_KEY = 'Baseline (prior)';
      // Separate regular models from the baseline entry so we can assign
      // deterministic colours and render the baseline last with a fixed style.
      const regularModels = Array.from(
        new Set([...Object.keys(roc), ...Object.keys(pr)])
      ).filter(m => m !== BASELINE_KEY);
      if (!regularModels.length && !pr[BASELINE_KEY] && !roc[BASELINE_KEY]) {
        el.innerHTML = '<p>No ROC/PR data available.</p>'; return;
      }

      const traces = [];
      regularModels.forEach((m, i) => {
        const c = palette[i % palette.length];
        const rocD = roc[m];
        const prD = pr[m];
        if (rocD && (rocD.fpr || []).length && (rocD.tpr || []).length) {
          traces.push({
            x: rocD.fpr, y: rocD.tpr, mode: 'lines', type: 'scatter',
            name: `${m} (AUC=${rocD.auc.toFixed(3)})`,
            legendgroup: m, showlegend: true,
            line: { color: c, width: 2 },
            xaxis: 'x', yaxis: 'y'
          });
        }
        if (prD && (prD.recall || []).length && (prD.precision || []).length) {
          traces.push({
            x: prD.recall, y: prD.precision, mode: 'lines', type: 'scatter',
            name: `${m} (AUPR=${prD.avg_precision.toFixed(3)})`,
            legendgroup: m, showlegend: true,
            line: { color: c, width: 2, dash: 'dash' },
            xaxis: 'x2', yaxis: 'y2'
          });
        }
      });

      // Baseline model (DummyClassifier prior) — grey dotted, own legend group.
      const bPrD = pr[BASELINE_KEY];
      const bRocD = roc[BASELINE_KEY];
      if (bRocD && (bRocD.fpr || []).length && (bRocD.tpr || []).length) {
        traces.push({
          x: bRocD.fpr, y: bRocD.tpr, mode: 'lines', type: 'scatter',
          name: `${BASELINE_KEY} (AUC=${bRocD.auc.toFixed(3)})`,
          legendgroup: BASELINE_KEY, showlegend: true,
          line: { color: '#888', width: 2, dash: 'dot' },
          xaxis: 'x', yaxis: 'y'
        });
      }
      if (bPrD && (bPrD.recall || []).length && (bPrD.precision || []).length) {
        traces.push({
          x: bPrD.recall, y: bPrD.precision, mode: 'lines', type: 'scatter',
          name: `${BASELINE_KEY} (AUPR=${bPrD.avg_precision.toFixed(3)})`,
          legendgroup: BASELINE_KEY, showlegend: true,
          line: { color: '#888', width: 2, dash: 'dashdot' },
          xaxis: 'x2', yaxis: 'y2'
        });
      }

      // ROC random diagonal.
      traces.push({
        x: [0, 1], y: [0, 1], mode: 'lines', type: 'scatter',
        name: 'Random (ROC)', line: { color: '#bbb', width: 1, dash: 'dot' },
        hoverinfo: 'skip', showlegend: false,
        xaxis: 'x', yaxis: 'y'
      });

      // PR random baseline: horizontal line at prevalence.
      if (prevalence != null) {
        traces.push({
          x: [0, 1], y: [prevalence, prevalence], mode: 'lines', type: 'scatter',
          name: `Random PR (prevalence=${prevalence.toFixed(3)})`,
          line: { color: '#bbb', width: 1, dash: 'dot' },
          hoverinfo: 'skip', showlegend: true,
          xaxis: 'x2', yaxis: 'y2'
        });
      }

        Plotly.newPlot(
            divId, traces,
            {
          margin: { l: 56, r: 20, t: 56, b: 90 },
          xaxis: { domain: [0.0, 0.46], title: 'False positive rate', range: [0, 1] },
          yaxis: { title: 'True positive rate', range: [0, 1] },
          xaxis2: { domain: [0.54, 1.0], title: 'Recall', range: [0, 1], anchor: 'y2' },
          yaxis2: { title: 'Precision', range: [0, 1], anchor: 'x2' },
          legend: { orientation: 'h', x: 0.5, y: -0.18, xanchor: 'center' },
          annotations: [
            { text: 'ROC curve (pooled across folds)', x: 0.23, y: 1.08, xref: 'paper', yref: 'paper', showarrow: false, font: { size: 14 } },
            { text: 'Precision-Recall curve (pooled across folds)', x: 0.77, y: 1.08, xref: 'paper', yref: 'paper', showarrow: false, font: { size: 14 } }
          ]
            },
            { responsive: true }
        );
    }

    // ── Confusion matrix heatmap ─────────────────────────────────────────────
    function plotConfusionMatrix(divId, cm, modelName, showColorbar) {
        // Row-normalize by true class (cm rows are [True 0, True 1]) so each
        // row sums to 1 -- classes are imbalanced enough that raw counts make
        // the minority true-class row unreadable (both in the printed number
        // and in the color scale, which raw counts would otherwise let the
        // majority row dominate). Absolute counts are kept in parentheses.
        const norm = cm.map((row) => {
            const rowSum = row[0] + row[1];
            return rowSum > 0 ? row.map((v) => v / rowSum) : row.map(() => 0);
        });
        const text = cm.map((row, i) =>
            row.map((v, j) => norm[i][j].toFixed(2) + ' (' + v + ')')
        );
        Plotly.newPlot(
            divId,
            [{
                z: norm,
                zmin: 0,
                zmax: 1,
                x: ['Pred 0', 'Pred 1'],
                y: ['True 0', 'True 1'],
                type: 'heatmap',
                // Explicit light-to-dark stops instead of the named 'Blues'
                // preset -- Plotly's built-in direction was rendering low
                // fractions dark and high fractions light, the opposite of
                // what a reader expects from a confusion-matrix heatmap.
                colorscale: [
                    [0, '#f7fbff'],
                    [0.25, '#c6dbef'],
                    [0.5, '#6baed6'],
                    [0.75, '#2171b5'],
                    [1, '#08306b']
                ],
                // Every panel shares the same fixed 0-1 scale, so only the
                // first panel in a grid needs to show the colorbar -- one
                // shared legend instead of one redundant copy per model.
                showscale: !!showColorbar,
                colorbar: { title: 'Fraction of true class', tickformat: '.0%' },
                text: text,
                texttemplate: '%{text}'
            }],
            {
                // Short title only -- the full "row-normalized by true
                // class" explanation lives once in the section header, not
                // repeated (and clipped) on every narrow grid panel.
                title: { text: modelName, font: { size: 13 } },
                margin: { l: 55, r: 15, t: 35, b: 40 },
                // Put True 0 on top so the diagonal (TN, TP) runs top-left
                // to bottom-right, matching standard confusion-matrix layout.
                yaxis: { autorange: 'reversed' }
            },
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
    run_datetime: str | None = None,
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
    html_generated_at = datetime.now(timezone.utc).isoformat()
    run_datetime_display = html.escape(run_datetime) if run_datetime else "N/A"
    html_generated_display = html.escape(html_generated_at)
    payload = [build_task_plot_payload(tr) for tr in task_results]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload_path = output_path.with_suffix(".json")
    payload_path.write_text(json.dumps(safe_json(payload)), encoding="utf-8")
    payload_file = payload_path.name
    model_card_html = _build_model_card_html(config, payload, run_datetime_display)

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
    .curves-row {{ display: grid; grid-template-columns: repeat(2, minmax(320px, 1fr)); gap: 12px; align-items: start; }}
    .curve-panel {{ min-width: 0; height: 400px; }}
    @media (max-width: 900px) {{
      .curves-row {{ grid-template-columns: 1fr; }}
    }}
    #spinner {{
      display: flex; flex-direction: column; align-items: center;
      justify-content: center; padding: 60px 0; gap: 16px; color: #666;
      font-size: 0.95rem;
    }}
    #spinner .spin {{
      width: 48px; height: 48px; border-radius: 50%;
      border: 5px solid #e0e0e0; border-top-color: #555;
      animation: spin 0.8s linear infinite;
    }}
    @keyframes spin {{ to {{ transform: rotate(360deg); }} }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <p style="color:#666;font-size:0.85rem;margin:4px 0 16px;">
    Results generated: <strong>{run_datetime_display}</strong> &nbsp;|&nbsp;
    Report rendered: <strong>{html_generated_display}</strong>
  </p>
{model_card_html}
  <div id="spinner"><div class="spin"></div><span>Loading report…</span></div>
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
      // Collected across every task card and rendered as one section at the
      // bottom of the report, instead of interleaved per-task-card.
      const allConfusionMatrices = [];
      for (const tr of reportData) {{
        const card = document.createElement('div');
        card.className = 'card';
        card.innerHTML = `<h2>Task: ${{tr.task}} (status=${{tr.status}})</h2><p>Primary metric: ${{tr.primary_metric}} | n_samples: ${{tr.n_samples}}</p>`;

        const p0 = document.createElement('div'); p0.id = `p0_${{tr.task}}`; p0.style.height = '420px';
        const p1 = document.createElement('div'); p1.id = `p1_${{tr.task}}`; p1.style.height = '320px';
        const p2 = document.createElement('div'); p2.id = `p2_${{tr.task}}`; p2.style.height = 'auto';
        const p3 = document.createElement('div'); p3.id = `p3_${{tr.task}}`; p3.style.height = '520px';
        const p4 = document.createElement('div'); p4.id = `p4_${{tr.task}}`; p4.style.height = '360px';
        card.appendChild(p0); card.appendChild(p1); card.appendChild(p2); card.appendChild(p3); card.appendChild(p4);

        const h3 = document.createElement('h3'); h3.textContent = 'Important parameters'; card.appendChild(h3);
        const ptab = document.createElement('div'); ptab.innerHTML = createParamsTable(tr.important_params); card.appendChild(ptab);
        container.appendChild(card);

        // ── Metric heatmap(s) ────────────────────────────────────────────────
        if (tr.task === 'regression') {{
          p3.style.display = 'none';
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

        // ── Model-wise fold distribution, repeated for every metric ──────────
        p1.innerHTML = '';
        p1.style.height = 'auto';
        const metricFoldData = tr.metric_fold_data || {{}};
        const metricKeys = Object.keys(metricFoldData);
        if (!metricKeys.length) {{
          p1.innerHTML = '<p>No fold-level metric data available.</p>';
        }} else {{
          const boxHeader = document.createElement('h3');
          boxHeader.textContent = 'Model-wise distribution per metric';
          p1.appendChild(boxHeader);
          const boxGrid = document.createElement('div');
          boxGrid.style.display = 'grid';
          boxGrid.style.gridTemplateColumns = 'repeat(auto-fit, minmax(340px, 1fr))';
          boxGrid.style.gap = '12px';
          p1.appendChild(boxGrid);
          // Create every container div up front so the CSS grid has already
          // settled on final column widths before any Plotly.newPlot call
          // measures its container -- rendering into a div mid-layout can
          // otherwise clip/mis-size the title on the first couple of plots.
          const boxDivByMetric = {{}};
          for (const metricName of metricKeys) {{
            const boxDiv = document.createElement('div');
            const safeMetricName = String(metricName).replace(/[^a-zA-Z0-9_]/g, '_');
            boxDiv.id = `box_${{tr.task}}_${{safeMetricName}}`;
            boxDiv.style.height = '320px';
            boxGrid.appendChild(boxDiv);
            boxDivByMetric[metricName] = boxDiv;
          }}
          for (const metricName of metricKeys) {{
            const md = metricFoldData[metricName] || {{ points: [], lines: [] }};
            const boxTraces = [];
            for (const modelName of tr.model_order) {{
              const yVals = (md.points || [])
                .filter(p => p.model_name === modelName)
                .map(p => p.score);
              boxTraces.push({{
                type: 'box', name: modelName, y: yVals,
                boxpoints: 'all', jitter: 0.35, pointpos: 0,
                marker: {{ size: 6, opacity: 0.75 }}, line: {{ width: 1 }}
              }});
            }}
            const foldLineTraces = (md.lines || [])
              .filter(fl => fl.x.length >= 2)
              .map(fl => ({{
                x: fl.x, y: fl.y, mode: 'lines+markers', type: 'scatter',
                name: `Fold ${{fl.outer_fold}}`, showlegend: false,
                marker: {{ size: 5, opacity: 0.65 }},
                line: {{ width: 1, color: '#666' }},
                hovertemplate: 'Fold ' + fl.outer_fold + '<br>%{{x}}: %{{y}}<extra></extra>'
              }}));
            const boxDiv = boxDivByMetric[metricName];
            Plotly.newPlot(
              boxDiv.id, [...boxTraces, ...foldLineTraces],
              {{
                title: {{ text: `Distribution of ${{metricName}}`, font: {{ size: 13 }} }},
                margin: {{ l: 55, r: 15, t: 45, b: 45 }},
                xaxis: {{ title: 'Model type' }},
                yaxis: {{ title: metricName }}
              }},
              {{ responsive: true }}
            );
          }}
        }}

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
          const thresholdPoints = [];
          for (let i = 0; i < thresholdRows.length; i++) {{
            const row = thresholdRows[i];
            const thr = Number(row.mean_threshold);
            if (!Number.isFinite(thr)) continue;
            const modelName = String(row.model_name || 'model');
            const c = colorByModel[modelName] || '#111';
            thresholdPoints.push({{
              modelName,
              color: c,
              rawThreshold: thr,
              displayThreshold: thr
            }});
          }}

          // Cluster nearby thresholds and apply small deterministic jitter so
          // per-model markers remain distinguishable without changing ordering.
          thresholdPoints.sort((a, b) => a.rawThreshold - b.rawThreshold);
          const CLUSTER_EPS = 0.012;
          const JITTER_STEP = 0.006;
          const MIN_X = 0.001;
          const MAX_X = 0.999;
          const clusters = [];
          for (const pt of thresholdPoints) {{
            const last = clusters[clusters.length - 1];
            if (!last || Math.abs(pt.rawThreshold - last[last.length - 1].rawThreshold) > CLUSTER_EPS) {{
              clusters.push([pt]);
            }} else {{
              last.push(pt);
            }}
          }}
          for (const cluster of clusters) {{
            const n = cluster.length;
            for (let i = 0; i < n; i++) {{
              const centered = i - ((n - 1) / 2);
              const jittered = cluster[i].rawThreshold + (centered * JITTER_STEP);
              cluster[i].displayThreshold = Math.max(MIN_X, Math.min(MAX_X, jittered));
            }}
          }}

          const MAX_LABEL_TIERS = 4;
          for (let i = 0; i < thresholdPoints.length; i++) {{
            const pt = thresholdPoints[i];
            const tier = i % MAX_LABEL_TIERS;
            const labelAy = -26 - (tier * 15);
            const xAnchor = (pt.displayThreshold > 0.9) ? 'right' : ((pt.displayThreshold < 0.1) ? 'left' : 'center');
            thresholdShapes.push({{
              type: 'line',
              x0: pt.displayThreshold,
              x1: pt.displayThreshold,
              y0: 0,
              y1: 1,
              yref: 'paper',
              line: {{ color: pt.color, width: 2, dash: 'dot' }}
            }});
            thresholdAnn.push({{
              x: pt.displayThreshold,
              y: 1,
              yref: 'paper',
              text: `${{pt.modelName}} thr=${{pt.rawThreshold.toFixed(3)}}`,
              showarrow: true,
              arrowhead: 2,
              arrowsize: 1,
              arrowwidth: 1,
              arrowcolor: pt.color,
              ax: 0,
              ay: labelAy,
              xanchor: xAnchor,
              yanchor: 'bottom',
              font: {{ size: 10, color: pt.color }}
            }});
          }}
          Plotly.newPlot(
            p2.id, probTraces,
            {{
              title: 'Predicted probability distribution by model (model thresholds shown)',
              xaxis: {{ title: 'P(class=1)' }}, yaxis: {{ title: 'Count' }},
              barmode: 'overlay',
              margin: {{ l: 60, r: 20, t: 95, b: 55 }},
              shapes: thresholdShapes,
              annotations: thresholdAnn
            }},
            {{ responsive: true }}
          );

          // Per-model confusion matrices are collected here and rendered
          // together in one section at the bottom of the report (see
          // allConfusionMatrices below the main render loop).
          const confByModel = tr.confusion_by_model || {{}};
          for (const modelName of Object.keys(confByModel)) {{
            allConfusionMatrices.push({{ task: tr.task, modelName, cm: confByModel[modelName] }});
          }}

          // ROC + PR curves with shared legend in one combined panel.
          plotRocPrCurves(p3.id, tr.roc_by_model || {{}}, tr.pr_by_model || {{}}, tr.pr_prevalence ?? null, palette);

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

        // ── Transformed (post-preprocessing) feature distributions ─────────
        const featureDist = tr.transformed_feature_distributions || {{}};
        const featureNames = Object.keys(featureDist);
        if (featureNames.length) {{
          const fdHeader = document.createElement('h3');
          fdHeader.textContent = 'Transformed feature distributions (after preprocessing)';
          card.appendChild(fdHeader);
          const fdNote = document.createElement('p');
          fdNote.style.color = '#666'; fdNote.style.fontSize = '0.85rem';
          fdNote.textContent = 'Sampled from the fitted preprocessor output (median-impute / log1p / percentile-clip as applicable, then StandardScaler); one-hot categorical columns omitted.';
          card.appendChild(fdNote);
          const fdGrid = document.createElement('div');
          fdGrid.style.display = 'grid';
          fdGrid.style.gridTemplateColumns = 'repeat(auto-fit, minmax(260px, 1fr))';
          fdGrid.style.gap = '10px';
          card.appendChild(fdGrid);
          // Two-pass render: create every container div before any
          // Plotly.newPlot call (same CSS-grid/Plotly timing fix used for
          // the per-metric box-plot grid above).
          const fdDivs = featureNames.map((name) => {{
            const safeName = String(name).replace(/[^a-zA-Z0-9_]/g, '_');
            const fdDiv = document.createElement('div');
            fdDiv.id = `fd_${{tr.task}}_${{safeName}}`;
            fdDiv.style.height = '220px';
            fdGrid.appendChild(fdDiv);
            return fdDiv;
          }});
          featureNames.forEach((name, i) => {{
            plotFeatureHistogram(fdDivs[i].id, featureDist[name], name);
          }});
        }}
      }}

      // ── Confusion matrices, gathered in one section at the report bottom ──
      if (allConfusionMatrices.length) {{
        const cmSection = document.createElement('div');
        cmSection.className = 'card';
        const cmSectionHeader = document.createElement('h2');
        cmSectionHeader.textContent = 'Confusion matrices (all models)';
        cmSection.appendChild(cmSectionHeader);
        const cmSectionNote = document.createElement('p');
        cmSectionNote.style.color = '#666';
        cmSectionNote.style.fontSize = '0.85rem';
        cmSectionNote.textContent = 'Row-normalized by true class -- each cell shows fraction of that row (absolute count in parentheses).';
        cmSection.appendChild(cmSectionNote);
        const cmGrid = document.createElement('div');
        cmGrid.style.display = 'grid';
        cmGrid.style.gridTemplateColumns = 'repeat(auto-fit, minmax(340px, 1fr))';
        cmGrid.style.gap = '12px';
        cmSection.appendChild(cmGrid);
        container.appendChild(cmSection);
        // Two-pass render: create every container div before any
        // Plotly.newPlot call (same CSS-grid/Plotly timing fix used for the
        // per-metric box-plot grid above).
        const cmDivs = allConfusionMatrices.map((entry) => {{
          const cmDiv = document.createElement('div');
          cmDiv.id = `cm_${{entry.task}}_${{entry.modelName}}`;
          cmDiv.style.height = '300px';
          cmGrid.appendChild(cmDiv);
          return cmDiv;
        }});
        allConfusionMatrices.forEach((entry, i) => {{
          plotConfusionMatrix(cmDivs[i].id, entry.cm, entry.modelName, i === 0);
        }});
      }}
    }}

    const spinner = document.getElementById('spinner');
    loadPayload(payloadFile)
      .then(renderReports)
      .catch((err) => {{
        const msg = document.createElement('div');
        msg.className = 'card';
        msg.innerHTML = `<h2>Report Load Error</h2><p>${{String(err)}}</p><p>If you opened this HTML via file://, serve the folder over HTTP so fetch() can read the JSON payload.</p>`;
        container.appendChild(msg);
      }})
      .finally(() => {{ spinner.style.display = 'none'; }});
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
    all_results: list[dict[str, Any]],
    out_dir: Path,
    verbose: bool = False,
    run_datetime: str | None = None,
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
        write_subset_html_report(
            key, task_results, report_path, run_datetime=run_datetime
        )
        vlog(verbose, f"Wrote HTML report: {report_path}", level="info")
