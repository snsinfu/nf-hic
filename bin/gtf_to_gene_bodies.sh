#!/usr/bin/env bash
#
# One gene-body interval per gene from a GTF, as sorted BED3 (chrom, start, end).
#
# Preference:
#   1. Rows with feature type "gene" (most faithful).
#   2. If none exist, collapse every attributed feature on (chrom, gene_id):
#      start = min(start)-1, end = max(end). This reconstructs gene bodies for
#      GTFs without "gene" rows (e.g. UCSC ncbiRefSeq/ensGene) and does not
#      over-count multi-isoform genes the way raw transcript rows would.
#
# "gene_id" is matched only at the start of a "; "-separated attribute token, so
# keys that merely end in gene_id (e.g. ref_gene_id) are not mistaken for it.
# Only uncompressed `*.gtf` input is accepted; other paths are rejected.

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $(basename "$0") <genes.gtf>" >&2
    exit 2
fi

gtf="$1"
case "$gtf" in
    *.gtf) ;;
    *)
        echo "ERROR: expected a .gtf file, got: $gtf" >&2
        exit 2
        ;;
esac

< "$gtf" awk -F'\t' -v OFS='\t' -v src="$gtf" '
    function get_gene_id(attrs,   i, n, parts, tok) {
        n = split(attrs, parts, ";")
        for (i = 1; i <= n; i++) {
            tok = parts[i]
            sub(/^[ \t]+/, "", tok)
            if (tok ~ /^gene_id[ \t]/) {
                sub(/^gene_id[ \t]+/, "", tok)
                gsub(/[";]/, "", tok)
                sub(/[ \t]+$/, "", tok)
                return tok
            }
        }
        return ""
    }
    /^#/ || NF < 9 { next }
    {
        c = $1; s = $4; e = $5
        id = get_gene_id($9)
        if (id != "") {
            k = c SUBSEP id
            if (!(k in cs) || s < cs[k]) cs[k] = s
            if (!(k in ce) || e > ce[k]) ce[k] = e
            if (!(k in cc)) { cc[k] = c; nall++ }
        }
        if ($3 == "gene") {
            gk = (id != "") ? c SUBSEP id : c SUBSEP s SUBSEP e
            if (!(gk in gs) || s < gs[gk]) gs[gk] = s
            if (!(gk in ge) || e > ge[gk]) ge[gk] = e
            if (!(gk in gc)) { gc[gk] = c; ngene++ }
        }
    }
    END {
        if (ngene > 0) {
            for (k in gs) print gc[k], gs[k] - 1, ge[k]
        } else {
            for (k in cs) print cc[k], cs[k] - 1, ce[k]
        }
        if (ngene == 0 && nall == 0) {
            print "ERROR: no gene bodies in " src ": need \"gene\" features or a gene_id attribute" > "/dev/stderr"
            exit 1
        }
    }
' | sort -k1,1 -k2,2n
