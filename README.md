# DNASeqAna

Small Python tool to detect periodic DNA repeats (motifs) from FASTA sequences. The tool parses FASTA input in chunks, builds position maps for dinucleotide tuples, detects periodic distance patterns and writes discovered repeats into a temporary SQLite database.

## Features
- Stream FASTA reading in chunks (`fasta_in_chunks`)
- Position indexing per dinucleotide (`process_sequence`)
- Periodicity detection from position lists (`calculate_difference`, `statistical_repeats`)
- Canonical motif determination (`canonical_dna_motif`)
- Temporary SQLite storage for results

## Requirements
- Python 3.8+
- Biopython

Install dependencies:

    pip install biopython

## Files
- `main.py` — main script with core logic and example run  
- `linkedList.py` — linked list implementation used by `main.py`  
- `test_data_split.fasta` — example/test FASTA used in `main.py`

## Usage
Run the script:

    python3 `main.py`

Default parameters inside `main.py`:
- `fasta_path` = `test_data_split.fasta`
- `chunk_size` = 5
- `motive_size` = 4
- `min_repeats` = 5
- `max_repeats` = 10

The script prints discovered motifs and stores results in a temporary SQLite database (created via `tempfile.NamedTemporaryFile`). An example `SELECT` at the end of `main.py` shows stored repeats.

## Implementation notes
- `fasta_in_chunks(fasta_path, chunk_size)` — yields lists of SeqRecord objects in chunks.  
- `process_sequence(seq)` — builds a dict of all dinucleotide keys (AA, AC, ... , TT) mapping to a `LinkedList` of positions.  
- `calculate_difference(ll)` — analyzes consecutive positions in a linked list to detect repeating distances (periods).  
- `statistical_repeats(base_tuple, seq, min_repeats, max_repeats)` — aggregates periodic findings and writes them to the SQLite table.  
- `canonical_dna_motif(seq)` — computes the lexicographically smallest rotation of a motif.  
- `reverse_complement(seq)` — returns the reverse complement (simple translate + reverse).

## TODOs & Known limitations
- Current implementation uses a custom `linkedList` module; consider replacing with built-in lists or optimized structures.  
- Chunking is fixed by sequence count; implement RAM-controlled chunking for large files.  
- Motif comparison / hashing and overlap handling need optimization.  
- `ProcessPoolExecutor` usage is present but commented out — adapt for parallel processing and large FASTA files.  
- Some function names contain minor typos (e.g., `lenght`) — review and fix unit tests.

## Example output
The script prints motif sequences and inserts rows into the `repeats` table. A sample `SELECT` is executed in `main.py`:

    SELECT * FROM repeats ORDER BY motif

## License
Unspecified. Add a `LICENSE` file if needed.
