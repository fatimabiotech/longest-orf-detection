# Longest ORF Detection and Sequence Analysis

A Python tool that analyzes DNA sequences from FASTA files. For each sequence it calculates base composition and GC content, finds all open reading frames (ORFs), identifies the longest ORF, translates it into a protein, and calculates the protein's molecular weight. Results are printed on screen and saved to a CSV file.

## Background
An open reading frame (ORF) is a stretch of DNA that starts with a start codon (`ATG`) and ends with a stop codon (`TAA`, `TAG` or `TGA`) in the same reading frame. Finding the longest ORF is a common first step in gene annotation and in checking which part of a sequence is likely to code for a protein.

## Features
- **Batch processing:** analyzes every `.fasta` file in the `data` folder in one run
- **Sequence statistics:** length, count of A, T, G and C, and GC content (%)
- **ORF detection:** finds every ORF from an `ATG` to the first in-frame stop codon
- **Longest ORF:** reports its coordinates, DNA length (bp) and protein length (amino acids)
- **Translation:** converts the longest ORF into its protein sequence
- **Molecular weight:** calculates the weight of the translated protein using Biopython
- **CSV export:** saves one row per file to `project_sequence_analysis.csv`

## Repository Structure
```
longest-orf-detection/
├── longest_orf_finding.py     # main script
├── data/
│   ├── PTEN_human_sequence.fasta
│   └── TP53_human_sequence.fasta
├── requirements.txt
└── README.md
```

## Requirements
- Python 3.x
- [Biopython](https://biopython.org/)

Install with:
```bash
pip install biopython
```

## How to Run
1. Clone or download this repository.
2. Put your FASTA files in the `data` folder (the human PTEN and TP53 sequences are included as examples).
3. Open a terminal in the repository folder and run:
```bash
python longest_orf_finding.py
```

## How It Works
1. The script reads each `.fasta` file in the `data` folder with Biopython (`SeqIO`).
2. `analyze_sequence()` counts A, T, G and C and calculates the GC content.
3. `find_orf_nucleotides()` searches for every `ATG` and moves forward in steps of three bases until it reaches a stop codon.
4. `find_longest_orf()` compares all ORFs, picks the longest one and translates it into a protein.
5. `find_protein_molecular_weight()` calculates the weight of that protein.
6. Everything is printed and written to `project_sequence_analysis.csv`.

## Output
For each file, the script prints the sequence ID, description, length, base counts, GC content, longest ORF coordinates, DNA length, protein length, protein sequence and molecular weight.

The CSV file contains these columns: `file_stem`, `id`, `gc_content`, `length`, `description`, `A`, `G`, `C`, `T`, `longest_coordinates`, `longest_dna`, `longest_dna_length`, `longest_protein`, `longest_protein_length`, `total_orfs`, `all_orfs`, `protein_MW`.

**Example output (replace with your real results):**
```
PTEN_human_sequence.fasta
  Coordinates: ___ to ___ (nucleotide positions)
  DNA length: ___ bp
  Protein length: ___ amino acids
```

## Limitations
- Only the forward strand is searched (the reverse complement is not analyzed)
- Only ORFs that end with a stop codon are reported
- Each FASTA file must contain exactly one sequence
- ORFs starting inside other ORFs are also counted, so `total_orfs` includes nested ORFs

## Future Improvements
- Search all six reading frames (add the reverse complement)
- Support FASTA files with multiple sequences
- Add a command-line option to choose the input folder and a minimum ORF length
- Plot ORF lengths with Matplotlib

## Author
Fatima Qasim, Biotechnology (MPhil) student interested in bioinformatics and computational biology.
GitHub: [fatimabiotech](https://github.com/fatimabiotech)
