#!/usr/bin/env python3
"""Export the Figure 3 matched, actual-threshold predictions for Frank."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def read_rows(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def export_rows(rows):
    selected = [r for r in rows if r["denominator"] == "figure3_matched"]
    keyed = {}
    for row in selected:
        key = (row["fromData"], row["toData"], float(row["effect_scale"]))
        if key in keyed or key[0] == key[1] or key[2] not in (1, 0.8):
            raise ValueError("Duplicate, self-comparison, or unexpected effect scale")
        keyed[key] = row
    studies = {key[0] for key in keyed}
    expected = {(a, b, scale) for a in studies for b in studies if a != b
                for scale in (1, 0.8)}
    if len(studies) != 10 or set(keyed) != expected:
        raise ValueError("Expected all 90 directed pairs at both effect scales")
    result = []
    # Preserve the original ordering by internal study IDs, before display-name mapping.
    pair_order = sorted({k[:2] for k in keyed}, key=lambda pair: (
        keyed[pair[0], pair[1], 1]["query"], keyed[pair[0], pair[1], 1]["target"]))
    for source, target in pair_order:
        full, atten = (keyed[source, target, scale] for scale in (1, 0.8))
        for column in ("target_n", "n_power_qtls", "n_observed_qtls", "observed",
                       "min_target_lod", "max_target_lod"):
            if full[column] != atten[column]:
                raise ValueError("Effect-scale rows disagree: " + column)
        n = int(full["n_power_qtls"])
        target_n = int(full["target_n"])
        observed = float(full["observed"])
        p_full = float(full["expected_actual_first_pass"])
        p_atten = float(atten["expected_actual_first_pass"])
        if n <= 0 or target_n <= 2 or n != int(full["n_observed_qtls"]):
            raise ValueError("Invalid or mismatched denominators")
        if not all(math.isfinite(p) and 0 <= p <= 1
                   for p in (observed, p_full, p_atten)) or p_atten > p_full + 1e-12:
            raise ValueError("Invalid probabilities")
        low, high = (float(full[c]) for c in ("min_target_lod", "max_target_lod"))
        if not (math.isfinite(low) and math.isfinite(high) and 0 < low <= high):
            raise ValueError("Invalid LOD range")
        count = round(n * observed)
        if abs(n * observed - count) > 1e-7:
            raise ValueError("Observed fraction does not recover an integer count")
        result.append(dict(
            source_study=source, target_study=target,
            target_sample_size=target_n, source_qtls_eligible=n,
            source_qtls_rediscovered=count,
            observed_rediscovery_fraction=full["observed"],
            predicted_rediscovery_fraction_full_effect=full["expected_actual_first_pass"],
            predicted_rediscovery_fraction_80pct_effect=atten["expected_actual_first_pass"],
            minimum_target_lod=full["min_target_lod"],
            maximum_target_lod=full["max_target_lod"]))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--reference", type=Path,
                        help="Check all original columns against a previous export")
    args = parser.parse_args()
    result = export_rows(read_rows(args.input))
    if args.reference:
        reference = read_rows(args.reference)
        prior = {(r["source_study"], r["target_study"]): r for r in reference}
        if len(reference) != 90 or len(prior) != 90 or set(prior) != {
                (r["source_study"], r["target_study"]) for r in result}:
            raise ValueError("Reference comparison keys differ")
        for row in result:
            for col, value in prior[row["source_study"], row["target_study"]].items():
                actual = row[col]
                if col in ("source_study", "target_study"):
                    equal = actual == value
                else:
                    equal = math.isclose(float(actual), float(value), rel_tol=1e-12, abs_tol=1e-14)
                if not equal:
                    raise ValueError("Reference mismatch: " + col)
    args.output_directory.mkdir(parents=True, exist_ok=False)
    path = args.output_directory / "rediscovery_predictions_comparisons.tsv"
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result[0]), delimiter="\t", lineterminator="\r\n")
        writer.writeheader()
        writer.writerows(result)
    metadata = {
        "input": str(args.input),
        "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "output_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "rows": len(result), "denominator": "figure3_matched",
        "required_generator": "iReviewFrankPower.R with same_direction_detection",
        "direction_note": "Direction is calculated upstream, not inferred or recalculated by this exporter.",
        "reference_verified": str(args.reference) if args.reference else None,
    }
    (args.output_directory / "rediscovery_predictions_provenance.json").write_text(
        json.dumps(metadata, indent=2) + "\n")
    print(path)


if __name__ == "__main__":
    main()
