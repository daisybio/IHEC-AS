#!/usr/bin/env python3
"""
Compute Pangolin reference splice site scores for SE and RI events.

Scores wild-type splice sites directly (no VCF / variant needed) using the
custom_usage.py API pattern from the Pangolin repo (Zeng & Li 2022, Genome Biology).
For each site: extract 10,001 bp centered on the junction, run through all 12
Pangolin ensemble models (3 members × 4 tissues), take P(splice site) per tissue.

Architecture notes (from pangolin/model.py):
  L=32, CL=2*sum(AR*(W-1))=10,000 → output length = input - 10,000 = 1 for SEQ_LEN=10,001
  i.e., for 10,001-bp input, only the central position is scored.
  Channels [1,4,7,10] = P(splice site) for Heart/Liver/Brain/Testis respectively.

Usage (from project root):
    mamba run -n ihec-as python scripts/compute_pangolin_scores.py [--batch-size 256]

Requires:
    see env.yml — pyfastx=2.3.0 (bioconda), pangolin pinned to commit 5cf94b8

Reads:
    processed_data/pangolin_events.csv         -- event_id, Event Type, ID (written by 03)
    data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    processed_data/5ss.fasta                  -- MaxEntScan 5'ss sequences (validation)
    processed_data/5ss_up.fasta               -- MaxEntScan upstream 5'ss sequences
    processed_data/3ss.fasta                  -- MaxEntScan 3'ss sequences

Outputs:
    processed_data/pangolin_scores.csv
"""

import argparse
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import importlib.resources
import numpy as np
import polars as pl
import torch
from tqdm import tqdm

try:
    import pyfastx
except ImportError:
    sys.exit("pyfastx not installed. Run: pip install pyfastx")

try:
    from pangolin.model import Pangolin, L, W, AR
except ImportError:
    sys.exit("pangolin not installed. Run: pip install pangolin (from tkzeng/Pangolin)")


# ── Constants ──────────────────────────────────────────────────────────────────

CONTEXT = 5000                  # bases each side of splice site
SEQ_LEN = 2 * CONTEXT + 1       # 10,001; model strips 5000 from each end → output = 1 position

# One-hot encoding: index 0=N, 1=A, 2=C, 3=G, 4=T
IN_MAP = np.asarray(
    [[0, 0, 0, 0],   # N
     [1, 0, 0, 0],   # A
     [0, 1, 0, 0],   # C
     [0, 0, 1, 0],   # G
     [0, 0, 0, 1]],  # T
    dtype=np.float32,
)

TISSUES = ["heart", "liver", "brain", "testis"]
# P(splice site) output channel per tissue (softmax second class from each output head)
TISSUE_CHANNELS = [1, 4, 7, 10]

# Per-filter invocation (PLAN §4.12e per-filter incremental outputs): the
# transcript_filter comes from the TRANSCRIPT_FILTER env var (Snakemake passes
# it per rule instance); interactive runs default to biotype_filtered. All
# event-set-specific inputs/outputs are suffixed with it.
TRANSCRIPT_FILTER = os.environ.get("TRANSCRIPT_FILTER", "biotype_filtered")

FASTA_PATH = Path("data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz")
OUTPUT_PATH = Path(f"processed_data/pangolin_scores_{TRANSCRIPT_FILTER}.csv")
RMATS_DIR   = Path(f"splicing_analysis/rmats/{TRANSCRIPT_FILTER}")

_BASE_LUT = np.zeros(256, dtype=np.int8)
for _char, _idx in (("A", 1), ("a", 1), ("C", 2), ("c", 2),
                     ("G", 3), ("g", 3), ("T", 4), ("t", 4)):
    _BASE_LUT[ord(_char)] = _idx


# ── Model loading ──────────────────────────────────────────────────────────────

def load_models(device: torch.device) -> list:
    """
    Load 12 Pangolin models: for i in [0,2,4,6] (tissue), j in [1,2,3] (ensemble).
    Weight path: pangolin/models/final.{j}.{i}.3.v2
    Returned list order: [tissue0_j1, tissue0_j2, tissue0_j3, tissue1_j1, ...]
    i.e. models[t*3 : (t+1)*3] = 3-member ensemble for tissue t.
    """
    models = []
    for i in [0, 2, 4, 6]:          # tissue index
        for j in range(1, 4):        # ensemble member
            m = Pangolin(L, W, AR)
            weight_path = str(
                importlib.resources.files("pangolin").joinpath(f"models/final.{j}.{i}.3.v2")
            )
            weights = torch.load(weight_path, map_location=device, weights_only=False)
            m.load_state_dict(weights)
            m.eval()
            m.to(device)
            models.append(m)
    return models


# ── Sequence encoding ──────────────────────────────────────────────────────────

def _encode(seq: str) -> np.ndarray:
    """DNA string → one-hot (4, L) float32. N → all-zero row."""
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    return IN_MAP[_BASE_LUT[arr]].T   # (4, L)


def encode_site(seq: str, strand: str) -> np.ndarray:
    """
    Encode a SEQ_LEN-bp string for Pangolin, applying RC for '-' strand.
    Pads with N if len(seq) < SEQ_LEN (e.g. near chromosome boundaries).
    Returns (4, SEQ_LEN) float32.
    Note: RC is applied to the input sequence (not the output), matching the
    custom_usage.py convention: (5 - seq_encoded) % 5 == RC.
    """
    if len(seq) < SEQ_LEN:
        pad_left  = (SEQ_LEN - len(seq)) // 2
        pad_right = SEQ_LEN - len(seq) - pad_left
        seq = "N" * pad_left + seq + "N" * pad_right
    enc = _encode(seq[:SEQ_LEN])    # (4, SEQ_LEN)
    if strand == "-":
        # Reverse complement: flip both axes (A↔T rows, C↔G rows, reverse time axis)
        enc = enc[::-1, ::-1].copy()
    return enc


# ── FASTA extraction ───────────────────────────────────────────────────────────

def extract_context(fa: "pyfastx.Fasta", chrom: str, pos: int) -> str:
    """
    Extract SEQ_LEN bp of + strand sequence centered on 0-based position pos.
    Pads with N where the window exceeds chromosome boundaries.
    Strand correction (RC) is done in encode_site, not here.
    """
    start = pos - CONTEXT
    end   = pos + CONTEXT + 1   # exclusive
    chrom_len = len(fa[chrom])

    clip_left  = max(0, -start)
    clip_right = max(0, end - chrom_len)
    seq = fa[chrom][max(0, start): min(chrom_len, end)].seq
    return "N" * clip_left + seq + "N" * clip_right


# ── Parallel batch encoding (worker pool) ───────────────────────────────────────
# extract_context + encode_site are pure-Python/CPU work; running them in worker
# processes keeps the GPU fed instead of sitting idle behind a serial Python loop.

_worker_fa = None

def _init_worker(fasta_path: Path) -> None:
    global _worker_fa
    _worker_fa = pyfastx.Fasta(str(fasta_path))

def _encode_one(site: dict) -> np.ndarray:
    seq = extract_context(_worker_fa, site["chrom"], site["pos"])
    return encode_site(seq, site["strand"])


# ── Coordinate parsing ─────────────────────────────────────────────────────────

def parse_suppa_sites(event_id: str, event_type: str):
    """
    Parse SUPPA event ID → (chrom, strand, sites_dict).
    sites_dict maps label → 0-based genomic position passed to Pangolin.

    Empirically validated positions (diag3.py across 200 events, GT/AG ≥ 0.97):

    SE  gene;SE:chr:e1-s2:e2-s3:strand
      + strand:
        upstream_donor      = e1      (GT at CONTEXT)
        alt_acceptor        = s2 - 1  (AG at CONTEXT-2:CONTEXT)
        alt_donor           = e2      (GT at CONTEXT)
        downstream_acceptor = s3 - 1  (AG at CONTEXT-2:CONTEXT)
      - strand (upstream/downstream swap in transcript order):
        upstream_donor      = s3 - 2  (GT at CONTEXT after RC)
        alt_acceptor        = e2 - 1  (AG at CONTEXT-2:CONTEXT after RC)
        alt_donor           = s2 - 2  (GT at CONTEXT after RC)
        downstream_acceptor = e1 - 1  (AG at CONTEXT-2:CONTEXT after RC)

    RI  gene;RI:chr:s1:e1-s2:e2:strand  (intron = [e1, s2))
      + strand: donor = e1,     acceptor = s2 - 1
      - strand: donor = s2 - 2, acceptor = e1 - 1

    These positions also align with R's 03-prepare-aggregation.Rmd MaxEntScan
    windows (no offset adjustment needed).
    """
    parts  = event_id.split(":")
    chrom  = parts[1]
    strand = parts[-1]

    if event_type == "SE":
        e1, s2 = map(int, parts[2].split("-"))
        e2, s3 = map(int, parts[3].split("-"))
        if strand == "+":
            sites = {
                "upstream_donor":      e1,
                "alt_acceptor":        s2 - 1,
                "alt_donor":           e2,
                "downstream_acceptor": s3 - 1,
            }
        else:
            sites = {
                "upstream_donor":      s3 - 2,
                "alt_acceptor":        e2 - 1,
                "alt_donor":           s2 - 2,
                "downstream_acceptor": e1 - 1,
            }
    elif event_type == "RI":
        # s1:e1-s2:e2
        s1 = int(parts[2])
        e1, s2 = map(int, parts[3].split("-"))
        e2 = int(parts[4])
        if strand == "+":
            sites = {
                "donor":    e1,
                "acceptor": s2 - 1,
            }
        else:
            sites = {
                "donor":    s2 - 2,
                "acceptor": e1 - 1,
            }
    else:
        return chrom, strand, {}

    return chrom, strand, sites


# ── Batched scoring ────────────────────────────────────────────────────────────

@torch.no_grad()
def score_batch(
    encodings: np.ndarray,
    models: list,
    device: torch.device,
) -> np.ndarray:
    """
    Score a batch of encoded sequences.

    encodings : (B, 4, SEQ_LEN)
    returns   : (B, 4) — mean P(splice site) across 3 ensemble members, per tissue.

    For SEQ_LEN=10,001 input, model output has shape (B, 12, 1); we read [:, channel, 0].
    """
    tensor = torch.tensor(encodings, dtype=torch.float32, device=device)
    tissue_scores = np.zeros((len(encodings), 4), dtype=np.float32)

    for t, channel in enumerate(TISSUE_CHANNELS):
        preds = []
        for model in models[t * 3: (t + 1) * 3]:
            out = model(tensor)                       # (B, 12, 1)
            preds.append(out[:, channel, 0].cpu().numpy())
        tissue_scores[:, t] = np.mean(preds, axis=0)

    return tissue_scores   # (B, 4)


# ── Validation: GT/AG dinucleotide check ──────────────────────────────────────

def check_canonical_dinucleotides(all_sites: list, fa: "pyfastx.Fasta") -> None:
    """
    After strand correction, donor sites should start with GT at CONTEXT:CONTEXT+2
    and acceptor sites should end with AG at CONTEXT-2:CONTEXT.
    Samples up to 2000 sites; warns if canonical fraction < 0.90.
    """
    import random
    donor_sites    = [s for s in all_sites if "donor"    in s["site_type"]]
    acceptor_sites = [s for s in all_sites if "acceptor" in s["site_type"]]

    def dinuc_frac(sites, expected, start, end):
        sample = random.sample(sites, min(2000, len(sites)))
        hits = 0
        for s in sample:
            seq = extract_context(fa, s["chrom"], s["pos"])
            enc = encode_site(seq, s["strand"])
            # decode central dinucleotide from one-hot
            window = enc[:, start:end]          # (4, 2)
            bases = "ACGT"
            dinuc = "".join(
                bases[np.argmax(window[:, i])] if window[:, i].sum() > 0 else "N"
                for i in range(end - start)
            )
            hits += int(dinuc == expected)
        return hits / len(sample)

    frac_gt = dinuc_frac(donor_sites,    "GT", CONTEXT,     CONTEXT + 2)
    frac_ag = dinuc_frac(acceptor_sites, "AG", CONTEXT - 2, CONTEXT)

    print(f"  Donor GT fraction:    {frac_gt:.3f} (expected ≥ 0.90)", flush=True)
    print(f"  Acceptor AG fraction: {frac_ag:.3f} (expected ≥ 0.90)", flush=True)

    if frac_gt < 0.90:
        warnings.warn(
            f"Donor GT fraction {frac_gt:.3f} < 0.90 — "
            "likely coordinate or strand error at donor sites"
        )
    if frac_ag < 0.90:
        warnings.warn(
            f"Acceptor AG fraction {frac_ag:.3f} < 0.90 — "
            "likely coordinate or strand error at acceptor sites"
        )


# ── Validation: MaxEntScan 9-mer / 23-mer comparison ─────────────────────────

def compare_to_maxentscan(all_sites: list, fa: "pyfastx.Fasta", id_map: dict) -> None:
    """
    Compare Pangolin-extracted splice site windows to MaxEntScan input sequences.

    MaxEntScan convention (from R export-sequences chunk):
      5'ss (donor):   9-mer  = genome[pos-3 : pos+6]  (3 exon + 6 intron)
      3'ss (acceptor): 23-mer = genome[pos-20: pos+3]  (20 intron + 3 exon)

    In Pangolin's SEQ_LEN-bp context (after strand correction):
      Splice site at index CONTEXT = 5000.
      9-mer  from Pangolin context: [CONTEXT-3 : CONTEXT+6]
      23-mer from Pangolin context: [CONTEXT-20: CONTEXT+3]

    MaxEntScan FASTA files:
      5ss.fasta     → SE alt_donor (e2/s2-2) + RI donor (e1/s2-2) keyed by ID
      5ss_up.fasta  → SE upstream_donor (e1/s3-2)              keyed by ID
      3ss.fasta     → SE alt_acceptor (s2-1/e2-1) + RI acceptor keyed by ID
    + strand positions; - strand positions after slash.

    id_map: event_id (str) → integer ID from pangolin_events.csv.
    Only events with a known ID can be cross-checked.
    """
    mes5_path    = Path(f"processed_data/5ss_{TRANSCRIPT_FILTER}.fasta")
    mes5up_path  = Path(f"processed_data/5ss_up_{TRANSCRIPT_FILTER}.fasta")
    mes3_path    = Path(f"processed_data/3ss_{TRANSCRIPT_FILTER}.fasta")

    if not (mes5_path.exists() and mes3_path.exists()):
        warnings.warn(
            "MaxEntScan FASTA files not found (processed_data/5ss.fasta, 3ss.fasta); "
            "run 03-prepare-aggregation.Rmd export-sequences chunk first"
        )
        return

    def load_mes_fasta(path):
        """Load MaxEntScan FASTA → dict[int_id] = seq_upper."""
        seqs = {}
        for name, seq in pyfastx.Fasta(str(path), build_index=False):
            try:
                seqs[int(name)] = seq.upper()
            except ValueError:
                pass
        return seqs

    mes5    = load_mes_fasta(mes5_path)
    mes5up  = load_mes_fasta(mes5up_path) if mes5up_path.exists() else {}
    mes3    = load_mes_fasta(mes3_path)

    # Which MaxEntScan file corresponds to each site_type?
    # (R code uses down5ss for SE alt_donor, up5ss for upstream_donor, event5ss for RI)
    mes_lookup = {
        # (event_type, site_type): mes_dict
        ("SE", "alt_donor"):           mes5,
        ("SE", "upstream_donor"):      mes5up,
        ("SE", "alt_acceptor"):        mes3,
        ("SE", "downstream_acceptor"): {},   # R uses 3ssdown for this; skip
        ("RI", "donor"):               mes5,
        ("RI", "acceptor"):            mes3,
    }

    mismatches, checked = 0, 0
    import random
    # Sample up to 5000 sites for speed
    sample = random.sample(all_sites, min(5000, len(all_sites)))

    for s in sample:
        ev_id = s["event_id"]
        ann_id = id_map.get(ev_id)
        if ann_id is None:
            continue
        mes_dict = mes_lookup.get((s["event_type"], s["site_type"]), {})
        if ann_id not in mes_dict:
            continue


        # Extract window from Pangolin context (after strand correction)
        seq = extract_context(fa, s["chrom"], s["pos"])
        enc = encode_site(seq, s["strand"])   # (4, SEQ_LEN)

        is_donor = "donor" in s["site_type"]
        if is_donor:
            ws, we = CONTEXT - 3, CONTEXT + 6    # 9-mer
        else:
            ws, we = CONTEXT - 20, CONTEXT + 3   # 23-mer

        # Decode only the needed window — decoding full SEQ_LEN is 10001 iterations
        window = enc[:, ws:we]  # (4, window_size)
        bases = "ACGT"
        pang_window = "".join(
            bases[np.argmax(window[:, i])] if window[:, i].sum() > 0 else "N"
            for i in range(we - ws)
        )

        if pang_window.upper() != mes_dict[ann_id]:
            mismatches += 1
        checked += 1

    if checked == 0:
        warnings.warn(
            "No events with matching IDs in both Pangolin coordinates and MaxEntScan files; "
            "cross-check skipped"
        )
        return

    mismatch_rate = mismatches / checked
    print(
        f"  MaxEntScan sequence cross-check: {checked} sites checked, "
        f"{mismatches} mismatches ({mismatch_rate:.2%})",
        flush=True,
    )
    if mismatch_rate > 0.05:
        raise ValueError(
            f"MaxEntScan vs Pangolin sequence mismatch rate {mismatch_rate:.1%} > 5%.\n"
            "Check: (1) FASTA path, (2) coordinate extraction off-by-one, "
            "(3) strand handling, (4) transcript_filter mismatch in ID lookup."
        )


# ── Event loading ──────────────────────────────────────────────────────────────

def load_events() -> tuple[list, dict]:
    """
    Load all SE + RI events.
    Returns:
      all_sites: list of site dicts (event_id, event_type, site_type, chrom, pos, strand)
      id_map:    event_id → integer ID (empty if pangolin_events.csv not found)
    """
    all_sites: list = []
    id_map:    dict = {}

    events_csv = Path(f"processed_data/pangolin_events_{TRANSCRIPT_FILTER}.csv")
    if events_csv.exists():
        df = pl.read_csv(events_csv)
        for row in df.iter_rows(named=True):
            event_id   = str(row["event_id"])
            event_type = str(row["Event Type"])
            id_map[event_id] = int(row["ID"])
            chrom, strand, sites = parse_suppa_sites(event_id, event_type)
            if not sites:
                continue
            for site_type, pos in sites.items():
                all_sites.append({
                    "event_id":   event_id,
                    "event_type": event_type,
                    "site_type":  site_type,
                    "chrom":      chrom,
                    "pos":        pos,
                    "strand":     strand,
                })
        return all_sites, id_map

    # Fallback: read from PSI file headers (no ID mapping, no MaxEntScan cross-check)
    warnings.warn(
        "processed_data/pangolin_events.csv not found; falling back to PSI headers "
        "(no integer ID mapping, MaxEntScan cross-check skipped)"
    )
    for event_type in ["SE", "RI"]:
        psi_file = RMATS_DIR / f"event_{event_type}.psi"
        if not psi_file.exists():
            warnings.warn(f"PSI file not found: {psi_file}")
            continue
        seen = set()
        with open(psi_file) as fh:
            for line in fh:
                event_id = line.split("\t")[0].strip()
                # Skip header / UUID lines that lack the expected semicolon
                if not event_id or ";" not in event_id or event_id in seen:
                    continue
                seen.add(event_id)
                chrom, strand, sites = parse_suppa_sites(event_id, event_type)
                if not sites:
                    continue
                for site_type, pos in sites.items():
                    all_sites.append({
                        "event_id":   event_id,
                        "event_type": event_type,
                        "site_type":  site_type,
                        "chrom":      chrom,
                        "pos":        pos,
                        "strand":     strand,
                    })

    return all_sites, id_map


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--skip-validation", action="store_true",
                        help="Skip GT/AG and MaxEntScan cross-checks (faster, for re-runs)")
    parser.add_argument("--n-events", type=int, default=None,
                        help="Score only the first N events (smoke test)")
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}", flush=True)
    if device.type == "cuda":
        print(f"  GPU: {torch.cuda.get_device_name(0)} "
              f"({torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB)", flush=True)

    print("Loading 12 Pangolin models...", flush=True)
    models = load_models(device)
    print("  models loaded.", flush=True)

    print(f"Opening FASTA: {FASTA_PATH}", flush=True)
    if not FASTA_PATH.exists():
        sys.exit(f"FASTA not found: {FASTA_PATH}")
    fa = pyfastx.Fasta(str(FASTA_PATH))
    print(f"  FASTA opened ({len(fa)} sequences).", flush=True)

    print("Loading event coordinates...", flush=True)
    all_sites, id_map = load_events()
    if args.n_events is not None:
        all_sites = all_sites[: args.n_events * 4]
        print(f"  --n-events {args.n_events}: truncated to {len(all_sites)} sites", flush=True)
    print(f"  {len(all_sites)} total sites to score "
          f"({len(id_map)} events with integer ID for cross-checks)", flush=True)

    # ── Validation checks ───────────────────────────────────────────────────────
    if not args.skip_validation:
        print("Checking GT/AG dinucleotides at splice site positions...", flush=True)
        check_canonical_dinucleotides(all_sites, fa)
        print("  GT/AG check passed.", flush=True)

        if id_map:
            print("Comparing Pangolin windows to MaxEntScan input sequences...", flush=True)
            compare_to_maxentscan(all_sites, fa, id_map)
            print("  MaxEntScan cross-check passed.", flush=True)
        else:
            print("  (MaxEntScan cross-check skipped: no integer ID map)", flush=True)

    # ── Batched scoring ─────────────────────────────────────────────────────────
    n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 4))
    print(f"Scoring {len(all_sites)} sites with batch_size={args.batch_size}, "
          f"encoding with {n_workers} worker processes...", flush=True)
    results = []
    n_batches = -(-len(all_sites) // args.batch_size)  # ceiling division

    pool = ProcessPoolExecutor(
        max_workers=n_workers, initializer=_init_worker, initargs=(FASTA_PATH,)
    )
    for b_start in tqdm(range(0, len(all_sites), args.batch_size), total=n_batches, file=sys.stdout):
        batch = all_sites[b_start: b_start + args.batch_size]

        encodings = np.stack(
            list(pool.map(_encode_one, batch, chunksize=16)),
            axis=0,
        )  # (B, 4, SEQ_LEN)

        ts = score_batch(encodings, models, device)   # (B, 4)

        for i, s in enumerate(batch):
            row_ts = ts[i]
            mean_p = float(np.mean(row_ts))
            results.append({
                "event_id":       s["event_id"],
                "event_type":     s["event_type"],
                "site_type":      s["site_type"],
                "tissue_heart":   float(row_ts[0]),
                "tissue_liver":   float(row_ts[1]),
                "tissue_brain":   float(row_ts[2]),
                "tissue_testis":  float(row_ts[3]),
                "mean_usage":     mean_p,
                "tissue_cv":      (float(np.std(row_ts) / mean_p)
                                   if mean_p > 0 else float("nan")),
            })

    pool.shutdown()

    # ── Write output ─────────────────────────────────────────────────────────────
    print("Writing output...", flush=True)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(results).write_csv(OUTPUT_PATH)
    print(f"Saved {len(results)} scores → {OUTPUT_PATH}", flush=True)


if __name__ == "__main__":
    main()
