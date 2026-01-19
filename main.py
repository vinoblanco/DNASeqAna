import itertools
import sqlite3, tempfile
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from linkedList import LinkedList
from typing import List, Iterator

try:
    import psutil
except ImportError:
    psutil = None

def _get_available_memory_bytes() -> int:
    """
    Return available RAM in bytes. Prefer psutil if installed,
    otherwise read /proc/meminfo (Linux).
    """
    if psutil:
        return int(psutil.virtual_memory().available)
    #Fallback for Linux
    meminfo_path = "/proc/meminfo"
    if os.path.exists(meminfo_path):
        with open(meminfo_path, "r") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    parts = line.split()
                    # value is in kB
                    return int(parts[1]) * 1024
    #Last resort: assume 256 MB available
    return 256 * 1024 * 1024

#todo motive abgleichen
#todo optimieren das nicht komplette Sequenz mehrfach analysiert

def fasta_in_chunks(fasta_path: str,
                    max_ram_mb: int = None,
                    ram_fraction: float = 0.5,
                    avg_bytes_per_base: float = 1.5) -> Iterator[List]:
    """
    Yield lists of SeqRecord objects sized by available RAM.
    Parameters:
    - fasta_path: path to FASTA file
    - max_ram_mb: optional explicit max RAM (MB) to use for a chunk
    - ram_fraction: fraction of available RAM to consume (0 < ram_fraction <= 1)
    - avg_bytes_per_base: estimated bytes used per base in memory (tweak as needed)
    """
    if max_ram_mb is not None:
        allowed_bytes = int(max_ram_mb * 1024 * 1024)
    else:
        avail = _get_available_memory_bytes()
        allowed_bytes = int(avail * float(ram_fraction))

    # safety floor
    allowed_bytes = max(allowed_bytes, 1 * 1024 * 1024)  # at least 1 MB

    chunk = []
    chunk_bytes = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        rec_bases = len(record.seq)
        rec_estimated_bytes = int(rec_bases * avg_bytes_per_base)

        # If a single record is larger than allowed_bytes, yield it alone
        if rec_estimated_bytes >= allowed_bytes:
            if chunk:
                yield chunk
                chunk = []
                chunk_bytes = 0
            yield [record]
            continue

        # If adding this record would exceed the allowed budget, yield current chunk
        if chunk and (chunk_bytes + rec_estimated_bytes) > allowed_bytes:
            yield chunk
            chunk = []
            chunk_bytes = 0

        chunk.append(record)
        chunk_bytes += rec_estimated_bytes

    if chunk:
        yield chunk

#fasta einlesen (test)
def fast_single_record(fasta_path):
    chunk = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        chunk.append(record)
        if len(chunk) == chunk_size:
            return chunk
    return None

def process_sequence(seq_str: str):
    """
    In dieser Methode wird ein Dictionary erstellt, welches mit Basentupeln als Keys und den vorkommen dieser als
    Value befüllt wird.
    :param seq_str: die Sequenz, für welche das Dic erstellt wird
    :return: das Dictionary
    """
    base_tuple = {''.join(kombi): None for kombi in itertools.product(['A','C','G','T'], repeat=2)}
    #print(base_tuple)

    for i in range(len(seq_str) - 1):
        key = str(seq_str[i:i+2])
        if base_tuple[key] is None:
            ll = LinkedList()
            base_tuple[key] = ll
            ll.append(i)
        else:
            base_tuple[key].append(i)

    return base_tuple

#test für hashmap(manuell)
def print_hashmap(base_tuple):
    for key, value in base_tuple.items():
        try:
            count = value.lenght()
            value.iterate()
        except:
            count = 0
        print(f"{key}: {count}")

def statistical_repeats(base_tuple, seq_str: str, seq_id, min_repeats, max_repeats):
    """
    Werdet das Dictionary aus und schreibt die gefundenen Repeats in eine Datenbank
    :param base_tuple: das Dictionary
    :param seq_str: die Sequenz des Dictionary
    :param min_repeats: mindest Anzahl der Vorkommen, damit es als Repeat angesehen wird
    :param max_repeats: maximale Anzahl der Vorkommen, damit es als Repeat angesehen wird
    """
    repeats = []
    for key, value in base_tuple.items():
        if value is None or value.lenght() < min_repeats:
            continue
        for counter, firststart in calculate_difference(value):
            for period, count in counter.items():
                if count >= min_repeats - 1:
                    start = firststart
                    end = start + period * count

                    if end <= len(seq_str):
                        motif = str(canonical_dna_motif(seq_str[start:start + period]))
                        #print(motif)
                        repeats.append((seq_id, motif, start, period, count+1))
    return repeats

#todo motive vergleichen (hash vergleichen)
def calculate_difference(ll):
    """
    Wertet die Liked List aus, welche alle vorkommen eines Basentupels enthalten. Berechnet ob der abstand zwischen
    zwei Vorkommen periodisch ist.
    :param ll: Linked List mit den Vorkommen
    :return: wann das Vorkommen startet (firststart) und wie oft es auftritt (count)
    """
    last = 0
    firststart = 0
    counter = {}
    for n, n1 in ll.pairwise():
         current = n1.data - n.data
         if current != last:
             yield counter, firststart
             counter = {}
             counter[current] = 1
             firststart = n.data
             last = current
         else:
             counter[current] += 1
             last = current
    if counter:
        yield counter, firststart

#gibt das reverse complement von einer Sequenz zurück
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def canonical_dna_motif(seq):
    """
    Sortiert eine Sequenz lexikographisch und gibt das kleinste Ergebnis zurück
    :param seq: Sequenz
    :return: kleinste Sortierung
    """
    candidates = [seq[i:] + seq[:i] for i in range(len(seq))]
    return min(candidates)

def worker_process_chunk(records, min_repeats, max_repeats):
    """records: list of (id, seq_str). Return list of repeat rows to insert."""
    rows = []
    for rec_id, seq_str in records:
        base_tuple = process_sequence(seq_str)
        rows.extend(statistical_repeats(base_tuple, seq_str, rec_id , min_repeats, max_repeats))
    return rows

if __name__ == "__main__":
    fasta_path = "test_data_split.fasta"
    chunk_size = 5
    motive_size = 4
    min_repeats = 5
    max_repeats = 10

    tmp = tempfile.NamedTemporaryFile(suffix=".db")
    conn = sqlite3.connect(tmp.name)

    cur = conn.cursor()
    cur.execute("""
    CREATE TABLE repeats (
        seq_number text NOT NULL,
        motif text NOT NULL,
        start integer NOT NULL,
        period integer NOT NULL,
        repeat integer NOT NULL,
        UNIQUE (seq_number, motif))
    """)

    cur.execute("CREATE INDEX idx_motif ON repeats (motif)")
    conn.commit()

    with ProcessPoolExecutor(max_workers=4) as executor:
        conn = sqlite3.connect(tmp.name)
        cur = conn.cursor()
        for chunk in fasta_in_chunks(fasta_path, ram_fraction=0.4):
            # convert SeqRecord -> lightweight tuples to ensure picklability
            lightweight = [(rec.id, str(rec.seq)) for rec in chunk]
            future = executor.submit(worker_process_chunk, lightweight, min_repeats, max_repeats)
            rows = future.result()  # blocking per-chunk; or collect futures and handle asynchronously
            if rows:
                cur.executemany("INSERT OR IGNORE INTO repeats VALUES (?,?,?,?,?)", rows)
                conn.commit()

    #test output todo: in .txt Datei schreiben
    cur.execute("""SELECT * FROM repeats WHERE motif = "ACCCCTCAGGGT" ORDER BY motif""")
    for row in cur.fetchall():
        print(row)