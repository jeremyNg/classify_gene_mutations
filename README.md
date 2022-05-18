Classification of gene mutations
===

Simple Python-based implementation for the classification of mutations into indel (frameshift/in-frame), nonsense and missense substitutions as a teaching aid for the section on gene mutations in the A-level H2 biology syllabus. As such, this is not as a substitute for other more established variant prediction tools such as SIFT and VEP.

Implementation details
===

The basic logic of this module is shown below.

Usage
===
This module requires Python > 3.2 (due to unittest requirements).
After cloning this repository, one can simply do the `python __main__.py -normal [normal sequence] -mutant [mutant sequence]`. The codon table is as follow:

codon_map = {"AAA":"A",
            "TTT":"C",
            "ATG":"M",
            "CTA":"G",
            "CTG":"G",
            "AAT":"W"}

All codons not present in the codon table is treated as a stop codon.

Likewise, the chemical property map is as follow:

amino_acid_map = {"A":"uncharged",
                   "G": "uncharged",
                   "C": "polar",
                   "M": "polar",
                   "W": "hydrophobic"}
