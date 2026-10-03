import gzip
import math
import os


def _open_vcf(filename):
    if os.fspath(filename).lower().endswith(".gz"):
        return gzip.open(filename, "rt", encoding="utf-8")
    return open(filename, "r", encoding="utf-8")


def analyze_vcf(filename):
    """Validate a VCF and calculate record, filter, allele, and genotype counts."""
    metrics = {
        "variant_records": 0,
        "snps": 0,
        "indels": 0,
        "other_variants": 0,
        "pass_records": 0,
        "filtered_records": 0,
        "missing_filter_records": 0,
        "samples": [],
        "sample_count": 0,
        "called_genotypes": 0,
        "missing_genotypes": 0,
        "heterozygous_genotypes": 0,
        "homozygous_reference_genotypes": 0,
        "homozygous_alternate_genotypes": 0,
        "mean_qual": None,
        "contig_count": 0,
    }
    errors = []
    contigs = set()
    qual_sum = 0.0
    qual_count = 0
    header_seen = False

    try:
        with _open_vcf(filename) as file_handle:
            for line_number, raw_line in enumerate(file_handle, start=1):
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith("##contig=<"):
                    contig_name = _metadata_value(line, "ID")
                    if contig_name:
                        contigs.add(contig_name)
                    continue
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    columns = line.split("\t")
                    if len(columns) < 8 or columns[:8] != [
                        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
                    ]:
                        errors.append(f"Invalid VCF column header on line {line_number}.")
                        break
                    if len(columns) > 8 and columns[8] != "FORMAT":
                        errors.append(f"Expected FORMAT column on line {line_number}.")
                        break
                    metrics["samples"] = columns[9:]
                    metrics["sample_count"] = len(columns[9:])
                    header_seen = True
                    continue
                if line.startswith("#"):
                    errors.append(f"Unexpected VCF header line on line {line_number}.")
                    break
                if not header_seen:
                    errors.append(f"VCF record appears before #CHROM header on line {line_number}.")
                    break

                fields = line.split("\t")
                expected_fields = 9 + metrics["sample_count"] if metrics["sample_count"] else 8
                if len(fields) != expected_fields:
                    errors.append(
                        f"VCF record on line {line_number} has {len(fields)} columns; "
                        f"expected {expected_fields}."
                    )
                    break
                if not fields[0] or not fields[3] or not fields[4]:
                    errors.append(f"VCF record on line {line_number} has an empty CHROM, REF, or ALT.")
                    break
                try:
                    position = int(fields[1])
                except ValueError:
                    errors.append(f"Invalid VCF position on line {line_number}: {fields[1]!r}.")
                    break
                if position < 1:
                    errors.append(f"VCF position on line {line_number} must be positive.")
                    break

                reference = fields[3].upper()
                alternates = fields[4].split(",")
                if any(not allele or allele == "." for allele in alternates):
                    errors.append(f"Invalid ALT allele on VCF line {line_number}.")
                    break
                metrics["variant_records"] += 1
                contigs.add(fields[0])

                if fields[6] == "PASS":
                    metrics["pass_records"] += 1
                elif fields[6] == ".":
                    metrics["missing_filter_records"] += 1
                else:
                    metrics["filtered_records"] += 1

                for alternate in alternates:
                    if len(reference) == 1 and len(alternate) == 1:
                        metrics["snps"] += 1
                    elif len(reference) != len(alternate):
                        metrics["indels"] += 1
                    else:
                        metrics["other_variants"] += 1

                if fields[5] != ".":
                    try:
                        quality = float(fields[5])
                    except ValueError:
                        errors.append(f"Invalid QUAL value on VCF line {line_number}.")
                        break
                    if not math.isfinite(quality):
                        errors.append(f"Non-finite QUAL value on VCF line {line_number}.")
                        break
                    qual_sum += quality
                    qual_count += 1

                if metrics["sample_count"]:
                    format_fields = fields[8].split(":")
                    genotype_index = format_fields.index("GT") if "GT" in format_fields else None
                    if genotype_index is None:
                        metrics["missing_genotypes"] += metrics["sample_count"]
                        continue
                    for sample_data in fields[9:]:
                        sample_fields = sample_data.split(":")
                        genotype = (
                            sample_fields[genotype_index]
                            if genotype_index < len(sample_fields)
                            else "."
                        )
                        alleles = genotype.replace("|", "/").split("/")
                        called_alleles = [allele for allele in alleles if allele != "."]
                        if not called_alleles or len(called_alleles) != len(alleles):
                            metrics["missing_genotypes"] += 1
                            continue
                        if any(not allele.isdigit() for allele in called_alleles):
                            errors.append(f"Invalid genotype on VCF line {line_number}.")
                            break
                        allele_indexes = [int(allele) for allele in called_alleles]
                        if any(index > len(alternates) for index in allele_indexes):
                            errors.append(f"Genotype allele index out of range on VCF line {line_number}.")
                            break
                        metrics["called_genotypes"] += 1
                        if len(set(allele_indexes)) > 1:
                            metrics["heterozygous_genotypes"] += 1
                        elif allele_indexes[0] == 0:
                            metrics["homozygous_reference_genotypes"] += 1
                        else:
                            metrics["homozygous_alternate_genotypes"] += 1
                    if errors:
                        break
    except (OSError, UnicodeError, EOFError) as exc:
        errors.append(f"Could not read VCF file: {exc}")

    if not header_seen and not errors:
        errors.append("VCF file is missing its #CHROM header.")
    metrics["contig_count"] = len(contigs)
    metrics["mean_qual"] = qual_sum / qual_count if qual_count else None
    validation = {
        "valid": not errors,
        "format": "vcf",
        "record_count": metrics["variant_records"],
        "errors": errors,
    }
    return {"validation": validation, "metrics": metrics}


def _metadata_value(line, field_name):
    prefix = f"{field_name}="
    for field in line.removeprefix("##contig=<").rstrip(">").split(","):
        if field.startswith(prefix):
            return field[len(prefix) :].strip('"')
    return None
