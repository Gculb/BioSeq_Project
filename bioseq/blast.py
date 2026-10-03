import math

from Bio.Blast import NCBIXML, NCBIWWW


SUPPORTED_PROGRAMS = {"blastn", "blastp", "blastx", "tblastn", "tblastx"}


def BLAST_search(
    sequence,
    program="blastn",
    database="nt",
    email=None,
    hitlist_size=50,
    expect=10.0,
):
    """Run a remote NCBI BLAST search and return one best-HSP summary per hit.

    The sequence may be a raw sequence or FASTA-formatted text. Query coverage
    measures the fraction of query positions spanned by the hit's best HSP.
    NCBI's public BLAST service applies usage limits and may take several
    minutes to return a result.
    """
    if not isinstance(sequence, str) or not sequence.strip():
        raise ValueError("sequence must be a non-empty sequence or FASTA string.")
    if program not in SUPPORTED_PROGRAMS:
        raise ValueError(
            f"Unsupported BLAST program {program!r}; expected one of "
            f"{', '.join(sorted(SUPPORTED_PROGRAMS))}."
        )
    if not isinstance(database, str) or not database.strip():
        raise ValueError("database must be a non-empty NCBI BLAST database name.")
    if not isinstance(hitlist_size, int) or isinstance(hitlist_size, bool) or hitlist_size < 1:
        raise ValueError("hitlist_size must be a positive integer.")
    if (
        not isinstance(expect, (int, float))
        or isinstance(expect, bool)
        or not math.isfinite(expect)
        or expect <= 0
    ):
        raise ValueError("expect must be greater than zero.")

    request_options = {
        "format_type": "XML",
        "hitlist_size": hitlist_size,
        "expect": expect,
    }
    if email is not None:
        email = str(email).strip()
        if not email:
            raise ValueError("email must be a non-empty address when supplied.")
        request_options["email"] = email

    response = NCBIWWW.qblast(
        program, database.strip(), sequence.strip(), **request_options
    )
    try:
        record = NCBIXML.read(response)
    finally:
        response.close()

    query_length = record.query_length
    hits = []
    for alignment in record.alignments:
        if not alignment.hsps:
            continue
        best_hsp = min(alignment.hsps, key=lambda hsp: hsp.expect)
        aligned_length = best_hsp.align_length
        query_span = max(0, best_hsp.query_end - best_hsp.query_start + 1)
        hits.append(
            {
                "hit_id": alignment.hit_id,
                "accession": alignment.accession,
                "description": alignment.hit_def,
                "percent_identity": (
                    best_hsp.identities / aligned_length * 100
                    if aligned_length
                    else 0
                ),
                "query_coverage": (
                    query_span / query_length * 100 if query_length else 0
                ),
                "evalue": best_hsp.expect,
                "bit_score": best_hsp.bits,
            }
        )

    return hits


blast_search = BLAST_search
