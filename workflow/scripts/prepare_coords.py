"""Validate a user-supplied coordinate sheet against the VCF sample list.

Reads a tab-delimited sheet with a header row and writes two files:

  coords      headerless TSV of sample<TAB>longitude<TAB>latitude, restricted to
              samples present in the VCF and in VCF order. This is the format
              every downstream R script expects, so the column order here is
              load-bearing: longitude first, then latitude.
  no_coords   samples present in the VCF but lacking usable coordinates, one per
              line. Consumed by the rule that builds the analysis VCF.

Accepted header spellings are deliberately generous (lon/long/x, lat/y,
sample/sample_id/...) because this is the one file every user has to write by
hand, but the values are checked strictly: a latitude of 120 is a mistake, not
data, and silently transposing it produces a map that looks plausible and is
wrong.
"""

# Send stdout and stderr to the rule's log file. Snakemake does this for R
# scripts (via sink) but not for python ones, which otherwise leaves the log
# empty and the traceback buried in the master log.
import sys as _sys
if "snakemake" in globals() and getattr(snakemake, "log", None):
    _sys.stdout = _sys.stderr = open(snakemake.log[0], "w")

import csv
import sys

SAMPLE_ALIASES = {
    "sample_id", "sample", "sample_name", "samples", "id", "indv", "individual", "name",
}
LON_ALIASES = {"longitude", "long", "lon", "lng", "x"}
LAT_ALIASES = {"latitude", "lat", "y"}

MISSING = {"", "na", "n/a", "nan", "none", "null", "."}


def die(message):
    sys.exit(f"ERROR (coordinate sheet): {message}")


def find_column(header, aliases, what, path):
    """Index of the first header field matching one of aliases."""
    normalized = [h.strip().lower().lstrip("*") for h in header]
    for i, name in enumerate(normalized):
        if name in aliases:
            return i
    die(
        f"no {what} column found in {path}\n"
        f"  header seen     : {header}\n"
        f"  accepted names  : {', '.join(sorted(aliases))}\n"
        "The sheet must be TAB-delimited and have a header row."
    )


def read_sheet(path):
    """{sample: (lon, lat)} plus the list of samples whose coordinates were blank."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))

    rows = [r for r in rows if r and any(field.strip() for field in r)]
    if not rows:
        die(f"{path} is empty")

    header = rows[0]
    if len(header) < 3:
        die(
            f"{path} looks like it is not tab-delimited: the header row split into "
            f"{len(header)} field(s): {header}"
        )

    i_sample = find_column(header, SAMPLE_ALIASES, "sample", path)
    i_lon = find_column(header, LON_ALIASES, "longitude", path)
    i_lat = find_column(header, LAT_ALIASES, "latitude", path)

    coords = {}
    blank = []
    for lineno, row in enumerate(rows[1:], start=2):
        if len(row) <= max(i_sample, i_lon, i_lat):
            die(f"{path} line {lineno}: expected at least {max(i_sample, i_lon, i_lat) + 1} columns, got {len(row)}")

        sample = row[i_sample].strip()
        if not sample:
            die(f"{path} line {lineno}: empty sample name")
        if sample in coords or sample in blank:
            die(f"{path} line {lineno}: sample {sample!r} appears more than once")

        raw_lon, raw_lat = row[i_lon].strip(), row[i_lat].strip()
        if raw_lon.lower() in MISSING or raw_lat.lower() in MISSING:
            blank.append(sample)
            continue

        try:
            lon, lat = float(raw_lon), float(raw_lat)
        except ValueError:
            die(
                f"{path} line {lineno}: sample {sample!r} has non-numeric coordinates "
                f"(longitude={raw_lon!r}, latitude={raw_lat!r})"
            )

        if lon == 0 and lat == 0:
            blank.append(sample)
            continue

        coords[sample] = (lon, lat)

    if not coords:
        die(f"{path} contains no usable coordinates")

    check_ranges(coords, path)
    return coords, blank


def check_ranges(coords, path):
    """Reject out-of-range values, and catch the classic lat/long transposition."""
    bad_lat = {s: v for s, v in coords.items() if abs(v[1]) > 90}
    bad_lon = {s: v for s, v in coords.items() if abs(v[0]) > 180}

    if bad_lat and all(abs(v[0]) <= 90 for v in coords.values()):
        example = next(iter(bad_lat.items()))
        die(
            f"{path}: the longitude and latitude columns look swapped.\n"
            f"  e.g. sample {example[0]!r} has longitude={example[1][0]}, latitude={example[1][1]}\n"
            "  latitude must be within [-90, 90]. Check the header names against the values."
        )

    for label, offenders, limit in (
        ("latitude", bad_lat, 90),
        ("longitude", bad_lon, 180),
    ):
        if offenders:
            shown = ", ".join(f"{s} ({v[1] if label == 'latitude' else v[0]})"
                              for s, v in list(offenders.items())[:5])
            die(
                f"{path}: {len(offenders)} sample(s) have {label} outside [-{limit}, {limit}]: "
                f"{shown}"
            )


def main():
    samples_file = snakemake.input["samples"]
    out_coords = snakemake.output["coords"]
    out_no_coords = snakemake.output["no_coords"]
    sheet = snakemake.params["coords_sheet"]
    use_coords = bool(snakemake.params["use_coords"])
    drop_missing = bool(snakemake.params["drop_samples_without_coords"])

    with open(samples_file) as fh:
        vcf_samples = [line.strip() for line in fh if line.strip()]

    # No coordinates at all: write empty files so the DAG shape stays the same.
    # Nothing downstream reads them, but the analysis-VCF rule still consumes
    # no_coords (empty => copy the VCF through unchanged).
    if not use_coords:
        open(out_coords, "w").close()
        open(out_no_coords, "w").close()
        print("use_coords is false — wrote empty coordinate files.")
        return

    coords, blank = read_sheet(sheet)

    matched = [s for s in vcf_samples if s in coords]
    missing = [s for s in vcf_samples if s not in coords]
    unused = sorted(set(coords) - set(vcf_samples))

    print(f"VCF samples                     : {len(vcf_samples)}")
    print(f"samples with coordinates        : {len(matched)}")
    print(f"samples without coordinates     : {len(missing)}")
    if blank:
        print(f"  ...of which blank in the sheet: {len(blank)}")
    if unused:
        print(
            f"sheet rows not in the VCF       : {len(unused)} (ignored) "
            f"e.g. {', '.join(unused[:5])}"
        )

    if not matched:
        die(
            f"none of the {len(vcf_samples)} VCF samples were found in {sheet}.\n"
            f"  first VCF samples  : {', '.join(vcf_samples[:5])}\n"
            f"  first sheet samples: {', '.join(sorted(coords)[:5])}\n"
            "Sample IDs must match the names in the VCF exactly."
        )

    if missing and not drop_missing:
        die(
            f"{len(missing)} VCF sample(s) have no coordinates and "
            "drop_samples_without_coords is false: "
            f"{', '.join(missing[:10])}{' ...' if len(missing) > 10 else ''}\n"
            "Add them to the sheet, or set drop_samples_without_coords: true to "
            "exclude them from the analysis."
        )

    with open(out_coords, "w") as fh:
        for sample in matched:
            lon, lat = coords[sample]
            fh.write(f"{sample}\t{lon}\t{lat}\n")

    with open(out_no_coords, "w") as fh:
        for sample in missing:
            fh.write(f"{sample}\n")


if __name__ == "__main__":
    main()
