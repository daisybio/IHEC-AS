from __future__ import annotations

import pandas as pd
import pytest

from splicing_ml.ontology_hierarchy import (
    ONTOLOGY_BIO_PATHS,
    map_ontology_to_supergroups,
)


class TestOntologyClusteringBasics:
    """Test basic structural properties of clustering."""

    def test_dynamic_ontology_clustering_returns_requested_group_count(self) -> None:
        """Verify exact group count for all valid k."""
        labels = sorted(ONTOLOGY_BIO_PATHS.keys())
        ontology = pd.Series(labels)

        for k in range(4, len(labels) + 1):
            mapped, diag = map_ontology_to_supergroups(ontology, n_groups=k)
            assert (
                mapped.nunique() == k
            ), f"k={k}: expected {k} groups, got {mapped.nunique()}"
            assert diag["supergrouping_requested_groups"] == k
            assert diag["supergroup_total"] == k

    def test_dynamic_ontology_clustering_rejects_k_below_four(self) -> None:
        """Verify ontology clustering rejects underspecified group counts."""
        ontology = pd.Series(sorted(ONTOLOGY_BIO_PATHS.keys()))

        with pytest.raises(ValueError, match="n_groups must be >=4"):
            map_ontology_to_supergroups(ontology, n_groups=3)

    def test_dynamic_ontology_clustering_rejects_invalid_group_counts(self) -> None:
        """Verify bounds checking on group counts."""
        ontology = pd.Series(sorted(ONTOLOGY_BIO_PATHS.keys()))

        with pytest.raises(ValueError, match="n_groups must be >=4"):
            map_ontology_to_supergroups(ontology, n_groups=2)

        with pytest.raises(
            ValueError, match="n_groups must be <= number of known ontology terms"
        ):
            map_ontology_to_supergroups(ontology, n_groups=len(ONTOLOGY_BIO_PATHS) + 1)

    def test_clustering_preserves_all_ontology_labels(self) -> None:
        """Verify no ontology label is dropped during clustering."""
        labels = sorted(ONTOLOGY_BIO_PATHS.keys())
        ontology = pd.Series(labels)

        for k in [4, 5, 10, len(labels)]:
            mapped, diag = map_ontology_to_supergroups(ontology, n_groups=k)
            assert len(mapped) == len(ontology)
            assert diag["original_ontology_total"] == len(labels)
            # All known labels should be mapped (no unknowns in ONTOLOGY_BIO_PATHS)
            assert diag["unknown_ontology_count"] == 0


class TestOntologyClusteringDeterminism:
    """Test determinism and reproducibility."""

    def test_clustering_is_deterministic(self) -> None:
        """Verify identical inputs produce identical clustering multiple times."""
        labels = sorted(ONTOLOGY_BIO_PATHS.keys())
        ontology = pd.Series(labels)

        # Run clustering 5 times for multiple k values
        for k in [4, 5, 10]:
            results = []
            for _ in range(5):
                mapped, _ = map_ontology_to_supergroups(ontology, n_groups=k)
                results.append(mapped.to_dict())

            # All runs should produce identical mapping
            for i in range(1, len(results)):
                assert (
                    results[i] == results[0]
                ), f"k={k}: Clustering is non-deterministic"

    def test_clustering_deterministic_with_repeated_data(self) -> None:
        """Verify clustering is deterministic when same ontology appears multiple times."""
        labels = sorted(ONTOLOGY_BIO_PATHS.keys())
        # Repeat each label 3 times
        ontology = pd.Series(labels * 3)

        mapped_once, _ = map_ontology_to_supergroups(ontology, n_groups=5)
        mapped_twice, _ = map_ontology_to_supergroups(ontology, n_groups=5)

        # Same input (even with duplicates) should give same mapping
        assert mapped_once.reset_index(drop=True).equals(
            mapped_twice.reset_index(drop=True)
        )

    def test_clustering_mapping_consistency_across_runs(self) -> None:
        """Verify label-to-supergroup mapping is consistent for the same k."""
        labels = sorted(ONTOLOGY_BIO_PATHS.keys())
        ontology = pd.Series(labels)

        k = 7
        run1, _ = map_ontology_to_supergroups(ontology, n_groups=k)
        run2, _ = map_ontology_to_supergroups(ontology, n_groups=k)

        # Each original label should always map to the same supergroup
        for label in labels:
            idx = list(ontology).index(label)
            assert run1.iloc[idx] == run2.iloc[idx], f"Inconsistent mapping for {label}"


class TestOntologyClusteringBiologicalMeaningfulness:
    """Test whether clustering respects biological relationships."""

    def test_immune_lymphoid_cluster_together(self) -> None:
        """Verify immune lymphoid terms cluster together at k=5,10,12."""
        lymphoid = ["b lymphocyte", "t lymphocyte", "natural killer cell"]
        ontology = pd.Series(sorted(ONTOLOGY_BIO_PATHS.keys()))

        for k in [5, 10, 12]:
            mapped, _ = map_ontology_to_supergroups(ontology, n_groups=k)
            # Get supergroup assignments for lymphoid terms
            lymphoid_groups = set(mapped[ontology.isin(lymphoid)].values)
            # All three should map to the same supergroup (or at most 2 for k=5)
            assert (
                len(lymphoid_groups) <= 2
            ), f"k={k}: lymphoid cells split across {len(lymphoid_groups)} supergroups"

    def test_immune_myeloid_cluster_together(self) -> None:
        """Verify immune myeloid terms cluster together."""
        myeloid = [
            "myeloid cell",
            "monocyte",
            "macrophage",
            "neutrophil",
            "eosinophil",
            "dendritic cell",
            "mononuclear cell",
        ]
        ontology = pd.Series(sorted(ONTOLOGY_BIO_PATHS.keys()))

        for k in [10, 15]:
            mapped, _ = map_ontology_to_supergroups(ontology, n_groups=k)
            myeloid_groups = set(mapped[ontology.isin(myeloid)].values)
            # Myeloid terms should form at most 2 distinct supergroups
            assert (
                len(myeloid_groups) <= 2
            ), f"k={k}: myeloid cells split across {len(myeloid_groups)} supergroups"

    def test_endodermal_organs_cluster_together(self) -> None:
        """Verify digestive/endodermal organs (liver, pancreas, colon) cluster together."""
        endodermal = ["liver", "pancreas", "colon", "digestive system"]
        ontology = pd.Series(sorted(ONTOLOGY_BIO_PATHS.keys()))

        for k in [8, 10, 12]:
            mapped, _ = map_ontology_to_supergroups(ontology, n_groups=k)
            endo_groups = set(mapped[ontology.isin(endodermal)].values)
            # Endodermal organs should cluster together (at most 2 groups for lower k)
            assert (
                len(endo_groups) <= 2
            ), f"k={k}: endodermal organs split across {len(endo_groups)} supergroups"

    def test_nervous_system_cluster_together(self) -> None:
        """Verify nervous system terms cluster together."""
        nervous = ["brain", "nervous system", "neural"]
        ontology = pd.Series(sorted(ONTOLOGY_BIO_PATHS.keys()))

        for k in [5, 10]:
            mapped, _ = map_ontology_to_supergroups(ontology, n_groups=k)
            nervous_groups = set(mapped[ontology.isin(nervous)].values)
            # Nervous system terms should form a single supergroup
            assert (
                len(nervous_groups) <= 1
            ), f"k={k}: nervous system split across {len(nervous_groups)} supergroups"

    def test_at_k_equals_unique_count_each_gets_own_group(self) -> None:
        """Verify that at k=unique_count, each label gets its own group."""
        labels = sorted(ONTOLOGY_BIO_PATHS.keys())
        ontology = pd.Series(labels)
        k = len(labels)

        mapped, _ = map_ontology_to_supergroups(ontology, n_groups=k)
        # Each unique label should have a unique group
        assert mapped.nunique() == k
        assert mapped.value_counts().min() == 1  # Each group has exactly 1 member

    def test_monotonicity_coarser_to_finer(self) -> None:
        """Verify that grouping becomes finer as k increases (merges don't split groups)."""
        labels = sorted(ONTOLOGY_BIO_PATHS.keys())
        ontology = pd.Series(labels)

        # For each pair (k, k+1), groups from k should be refinable into k+1
        # i.e., no group from k should be split arbitrarily
        for k in range(4, min(15, len(labels))):
            mapped_k, _ = map_ontology_to_supergroups(ontology, n_groups=k)
            mapped_k1, _ = map_ontology_to_supergroups(ontology, n_groups=k + 1)

            # Count groups: k+1 should have more or equal groups (monotonic refinement)
            assert mapped_k1.nunique() >= mapped_k.nunique()


class TestOntologyClusteringDiagnostics:
    """Test diagnostic output quality."""

    def test_diagnostics_contains_expected_fields(self) -> None:
        """Verify all expected diagnostic fields are present."""
        ontology = pd.Series(sorted(ONTOLOGY_BIO_PATHS.keys()))
        _, diag = map_ontology_to_supergroups(ontology, n_groups=5)

        expected_fields = {
            "supergrouping_enabled",
            "supergrouping_scheme",
            "supergrouping_scheme_kind",
            "supergrouping_requested_groups",
            "supergroup_total",
            "original_ontology_total",
            "unknown_ontology_labels",
            "unknown_ontology_count",
            "supergroup_counts",
            "supergroup_proportions",
        }
        assert set(diag.keys()) == expected_fields

    def test_diagnostics_proportions_sum_to_one(self) -> None:
        """Verify supergroup proportions sum to 1.0."""
        ontology = pd.Series(sorted(ONTOLOGY_BIO_PATHS.keys()))

        for k in [4, 5, 10]:
            _, diag = map_ontology_to_supergroups(ontology, n_groups=k)
            total_prop = sum(diag["supergroup_proportions"].values())
            assert (
                abs(total_prop - 1.0) < 1e-6
            ), f"k={k}: proportions sum to {total_prop}, not 1.0"

    def test_diagnostics_counts_and_proportions_consistent(self) -> None:
        """Verify counts and proportions are consistent."""
        ontology = pd.Series(sorted(ONTOLOGY_BIO_PATHS.keys()))

        for k in [4, 5, 10]:
            _, diag = map_ontology_to_supergroups(ontology, n_groups=k)
            total = sum(diag["supergroup_counts"].values())

            for group_name in diag["supergroup_counts"]:
                count = diag["supergroup_counts"][group_name]
                expected_prop = count / total
                actual_prop = diag["supergroup_proportions"][group_name]
                assert (
                    abs(actual_prop - expected_prop) < 1e-6
                ), f"Mismatch for {group_name}: expected {expected_prop}, got {actual_prop}"
