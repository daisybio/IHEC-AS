from __future__ import annotations

"""Ontology-to-supergroup mapping utilities for grouped CV.

Uses biology-guided hierarchical clustering to dynamically partition ontology
labels into k supergroupings for any k in [4, unique_ontology_count]. Clustering
is based on hierarchical biological paths that encode coarse-to-fine biological
proximity, ensuring stable and biologically meaningful splits across runs.
"""

from collections import Counter
from typing import Any

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.spatial.distance import squareform


__all__ = [
    "map_ontology_to_supergroups",
]


# Normalized ontology label -> hierarchical biological path.
# Paths are coarse-to-fine and encode biological proximity for distance-based
# clustering. The final leaf is the original ontology term.
ONTOLOGY_BIO_PATHS: dict[str, tuple[str, ...]] = {
    # Immune (lymphoid)
    "b lymphocyte": ("immune", "lymphoid", "b_lymphocyte"),
    "t lymphocyte": ("immune", "lymphoid", "t_lymphocyte"),
    "natural killer cell": ("immune", "lymphoid", "natural_killer"),
    # Immune (myeloid)
    "myeloid cell": ("immune", "myeloid", "myeloid_cell"),
    "monocyte": ("immune", "myeloid", "monocyte"),
    "macrophage": ("immune", "myeloid", "macrophage"),
    "neutrophil": ("immune", "myeloid", "neutrophil"),
    "eosinophil": ("immune", "myeloid", "eosinophil"),
    "dendritic cell": ("immune", "myeloid", "dendritic"),
    "mononuclear cell": ("immune", "myeloid", "mononuclear"),
    # Blood / hematopoietic lineage
    "hematopoietic cell": ("immune", "hematopoietic", "hematopoietic"),
    "erythroid lineage cell": ("immune", "hematopoietic", "erythroid"),
    "peripheral blood": ("immune", "hematopoietic", "peripheral_blood"),
    # Nervous
    "brain": ("nervous_system", "brain", "brain"),
    "nervous system": ("nervous_system", "general", "nervous_system"),
    "neural": ("nervous_system", "general", "neural"),
    # Digestive / hepatopancreatic / endodermal
    "digestive system": ("endodermal", "digestive", "digestive_system"),
    "colon": ("endodermal", "digestive", "colon"),
    "liver": ("endodermal", "hepatopancreatic", "liver"),
    "pancreas": ("endodermal", "hepatopancreatic", "pancreas"),
    "endoderm-derived structure": ("endodermal", "general", "endodermal_structure"),
    # Epithelial and mucosal
    "epithelial": ("epithelial_mucosal", "epithelial", "epithelial"),
    "mucosa": ("epithelial_mucosal", "mucosal", "mucosa"),
    # Mesoderm / connective / muscle
    "connective tissue cell": ("mesodermal", "connective", "connective_tissue"),
    "mesoderm-derived structure": ("mesodermal", "general", "mesodermal_structure"),
    "muscle": ("mesodermal", "muscle", "muscle"),
    # Organ-specific
    "kidney": ("organ_specific", "renal", "kidney"),
    "lung": ("organ_specific", "respiratory", "lung"),
    # Developmental / stem / extraembryonic
    "stem cell": ("developmental", "stem", "stem_cell"),
    "embryonic cell (metazoa)": ("developmental", "embryonic", "embryonic_cell"),
    "extraembryonic cell": ("developmental", "extraembryonic", "extraembryonic_cell"),
    "placenta": ("developmental", "extraembryonic", "placenta"),
    # Specialized / transformed
    "cancer cell line": ("specialized_transformed", "transformed", "cancer_cell_line"),
    "melanocyte": ("specialized_transformed", "specialized", "melanocyte"),
}


def _normalize_ontology_label(value: Any) -> str:
    """Internal helper for normalize ontology label."""
    return str(value).strip().lower()


def _common_prefix_length(a: tuple[str, ...], b: tuple[str, ...]) -> int:
    """Internal helper for common prefix length."""
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def _path_distance(a: tuple[str, ...], b: tuple[str, ...]) -> float:
    """Internal helper for path distance."""
    if a == b:
        return 0.0
    lcp = _common_prefix_length(a, b)
    # Tree-edit distance normalized to [0, 1].
    raw = (len(a) - lcp) + (len(b) - lcp)
    denom = max(1, len(a) + len(b))
    return float(raw / denom)


def _cluster_name_from_paths(paths: list[tuple[str, ...]], cluster_id: int) -> str:
    """Internal helper for cluster name from paths."""
    if not paths:
        return f"hc_cluster_{cluster_id}"
    prefix = list(paths[0])
    for p in paths[1:]:
        upto = _common_prefix_length(tuple(prefix), p)
        prefix = prefix[:upto]
        if not prefix:
            break
    if prefix:
        return "hc_" + "_".join(prefix)
    return f"hc_cluster_{cluster_id}"


def _build_hierarchical_mapping(n_groups: int) -> dict[str, str]:
    """Internal helper for build hierarchical mapping."""
    labels = sorted(ONTOLOGY_BIO_PATHS.keys())
    if n_groups < 4:
        raise ValueError(f"n_groups must be >=4, got {n_groups}")
    if n_groups > len(labels):
        raise ValueError(
            f"n_groups must be <= number of known ontology terms ({len(labels)}), got {n_groups}"
        )

    n = len(labels)
    dist = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = _path_distance(
                ONTOLOGY_BIO_PATHS[labels[i]], ONTOLOGY_BIO_PATHS[labels[j]]
            )
            dist[i, j] = d
            dist[j, i] = d

    condensed = squareform(dist, checks=False)
    z = linkage(condensed, method="average")
    # cut_tree gives exact n_clusters for valid bounds; fcluster(maxclust)
    # may return fewer than requested depending on distance plateaus.
    cluster_ids = cut_tree(z, n_clusters=n_groups).reshape(-1) + 1

    members: dict[int, list[str]] = {}
    for label, cid in zip(labels, cluster_ids):
        members.setdefault(int(cid), []).append(label)

    id_to_name: dict[int, str] = {}
    used_names: dict[str, int] = {}
    for cid, labs in sorted(members.items()):
        paths = [ONTOLOGY_BIO_PATHS[l] for l in labs]
        base_name = _cluster_name_from_paths(paths, cluster_id=cid)
        if base_name in used_names:
            used_names[base_name] += 1
            name = f"{base_name}_{used_names[base_name]}"
        else:
            used_names[base_name] = 1
            name = base_name
        id_to_name[cid] = name

    return {label: id_to_name[int(cid)] for label, cid in zip(labels, cluster_ids)}


def map_ontology_to_supergroups(
    ontology: pd.Series,
    n_groups: int,
) -> tuple[pd.Series, dict[str, Any]]:
    """Map ontology labels to biology-guided hierarchical supergroups.

    The number of output groups is dynamic and must be between 4 and the number
    of known ontology labels. Unknown labels are mapped to ``other``.
    """
    mapping = _build_hierarchical_mapping(n_groups)

    normalized = ontology.astype(str).map(_normalize_ontology_label)
    mapped = normalized.map(lambda x: mapping.get(x, "other"))

    original_counts = Counter(normalized.tolist())
    mapped_counts = Counter(mapped.tolist())
    total = max(1, int(len(mapped)))

    unknown_labels = sorted(
        [label for label in original_counts if label not in mapping]
    )

    diagnostics: dict[str, Any] = {
        "supergrouping_enabled": True,
        "supergrouping_scheme": f"hc{n_groups}",
        "supergrouping_scheme_kind": "hierarchical_dynamic",
        "supergrouping_requested_groups": int(n_groups),
        "supergroup_total": int(len(set(mapped.tolist()))),
        "original_ontology_total": int(len(original_counts)),
        "unknown_ontology_labels": unknown_labels,
        "unknown_ontology_count": int(
            sum(original_counts[label] for label in unknown_labels)
        ),
        "supergroup_counts": {
            key: int(val)
            for key, val in sorted(
                mapped_counts.items(), key=lambda kv: (-kv[1], kv[0])
            )
        },
        "supergroup_proportions": {
            key: float(val / total)
            for key, val in sorted(
                mapped_counts.items(), key=lambda kv: (-kv[1], kv[0])
            )
        },
    }

    return mapped.rename("ontology_supergroup"), diagnostics
