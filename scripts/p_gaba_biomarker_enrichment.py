#!/usr/bin/env python
"""Backward-compatible entry point for P GABA biomarker enrichment.

Delegates to :mod:`celltype_biomarker_enrichment` with P GABA defaults when
``--celltype_value`` and ``--gene_set_file`` are not supplied on the CLI.
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from celltype_biomarker_enrichment import main as _main

_DEFAULT_GENE_SET = _REPO_ROOT / "reference" / "p_gaba_genes_allen_yao_et_al.txt"


def main() -> None:
    """Inject P GABA defaults, then run :func:`celltype_biomarker_enrichment.main`.

    Returns:
        None
    """
    argv = sys.argv[1:]
    if not any(arg == "--celltype_value" or arg.startswith("--celltype_value=") for arg in argv):
        argv = ["--celltype_value", "P GABA", *argv]
    if not any(arg == "--gene_set_file" or arg.startswith("--gene_set_file=") for arg in argv):
        argv = ["--gene_set_file", str(_DEFAULT_GENE_SET), *argv]
    sys.argv = [sys.argv[0], *argv]
    _main()


if __name__ == "__main__":
    main()
