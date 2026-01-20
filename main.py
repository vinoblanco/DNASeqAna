import itertools
import sqlite3, tempfile
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from linkedList import LinkedList
from typing import List, Iterator, Union

try:
    import psutil
except ImportError:
    psutil = None

def _get_available_memory_bytes() -> int:
    """
    Gibt den verfügbaren RAM in Bytes zurück. Bevorzuge psutil, wenn installiert,
    ansonsten lese /proc/meminfo (Linux).
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
    Gibt Listen von SeqRecord Objekten zurück, die nach verfügbarem RAM dimensioniert sind.
    - fasta_path: Pfad zur FASTA Datei
    - max_ram_mb: maximaler RAM in MB (überschreibt ram_fraction, wenn gesetzt)
    - ram_fraction: Anteil des verfügbaren RAMs, der genutzt werden soll (0.0 - 1.0)
    - avg_bytes_per_base: geschätzte durchschnittliche Bytes pro Base in SeqRecord
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

def process_sequence(seq_str: str):
    """
    In dieser Methode wird ein Dictionary erstellt, welches mit Basentupeln als Keys und den Vorkommen dieser als
    Value befüllt wird.
    :param seq_str: die Sequenz, für welche das Dict erstellt wird
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

def statistical_repeats(base_tuple, seq_str: str, seq_id, min_repeats, max_repeats, motive_size):
    """
    Wertet das Dictionary aus und schreibt die gefundenen Repeats in eine Datenbank
    :param base_tuple: das Dictionary
    :param seq_str: die Sequenz des Dictionary als String
    :param seq_id: die ID der Sequenz
    :param min_repeats: mindest Anzahl der Vorkommen, damit es als Repeat angesehen wird
    :param max_repeats: maximale Anzahl der Vorkommen, damit es als Repeat angesehen wird
    :param motive_size: mindest Motivgröße
    :return: ein Array mit den gefundenen Repeats und deren Eigenschaften (Seq_ID, Motiv, Start, Periode, Anzahl)
    """
    repeats = []
    for key, value in base_tuple.items():
        if value is None or value.lenght() < min_repeats:
            continue
        for counter, first_start in calculate_difference(value, motive_size):
            for period, count in counter.items():
                if count >= min_repeats - 1:
                    start = first_start
                    end = start + period * count

                    if end <= len(seq_str):
                        motif = str(canonical_dna_motif(seq_str[start:start + period]))
                        #print(motif)
                        repeats.append((seq_id, motif, start, period, count+1))
    return repeats

#todo motive vergleichen (hash vergleichen)
def calculate_difference(ll, motive_size):
    """
    Wertet die Linked List aus, welche alle vorkommen eines Basentupels enthalten. Berechnet ob der Abstand zwischen
    zwei Vorkommen periodisch ist.
    :param ll: Linked List mit den Vorkommen des Basentupels
    :param motive_size: mindest Motivgröße
    :return: wann das Vorkommen startet (start) und wie oft es auftritt (count)
    """
    counter = {}
    start = None
    last = None

    for n, n1 in ll.pairwise():
        current = n1.data - n.data

        if last is None:
            last = current
            counter = {current: 1}
            start = n.data
            continue

        if current != last or current < motive_size:
            if counter:
                yield counter, start
            counter = {current: 1}
            start = n.data
            last = current
        else:
            counter[current] = counter.get(current, 0) + 1
            last = current

    if counter:
        yield counter, start

def reverse_complement(seq):
    """
    Gibt das reverse complement einer DNA Sequenz zurück
    :param seq: DNA Sequenz
    :return: reverse complement der DNA Sequenz
    """
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def canonical_dna_motif(seq):
    """
    Sortiert eine Sequenz lexikographisch und gibt das kleinste Ergebnis zurück
    :param seq: Sequenz
    :return: kleinste Sortierung
    """
    candidates = [seq[i:] + seq[:i] for i in range(len(seq))]
    return min(candidates)

def worker_process_chunk(records, min_repeats, max_repeats, motive_size):
    """
    Vorbereitet die Daten für die Verarbeitung in einem Prozess
    :param records: Sequenz Datensätze nur mit id und Sequenz String
    :param min_repeats: minimale Anzahl der Wiederholungen
    :param max_repeats: maximale Anzahl der Wiederholungen
    :param motive_size: mindest Motivgröße
    :return: Daten
    """
    rows = []
    for rec_id, seq_str in records:
        base_tuple = process_sequence(seq_str)
        repeats = statistical_repeats(base_tuple, seq_str, rec_id , min_repeats, max_repeats, motive_size)
        rows.extend(repeats)
    return rows

def write_repeats_to_txt(db: Union[str, sqlite3.Connection], output_path: str = "output.txt") -> None:
    """
    Schreibt die gefundenen Repeats aus der Datenbank in eine Textdatei.
    :param output_path: Pfad zur Ausgabedatei.
    :param db: Datenbankverbindung oder Pfad zur Datenbankdatei.
    """
    close_conn = False
    if isinstance(db, str):
        conn = sqlite3.connect(db)
        close_conn = True
    else:
        conn = db

    cur = conn.cursor()
    cur.execute("SELECT seq_number, motif, start, period, repeat FROM repeats ORDER BY motif ASC")
    rows = cur.fetchall()

    with open(output_path, "w", encoding="utf-8") as f:
        # header
        f.write("seq_number\tmotif\tstart\tperiod\trepeat\n")
        for row in rows:
            f.write("\t".join(str(col) for col in row) + "\n")

    if close_conn:
        conn.close()

if __name__ == "__main__":
    fasta_path = "test_data_split.fasta"
    motive_size = 4
    min_repeats = 3
    max_repeats = 10 #todo nutzen

    print("Starting processing...")

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
        futures = []
        print("Processing FASTA in chunks...")
        for chunk in fasta_in_chunks(fasta_path, ram_fraction=0.4):
            lightweight = [(rec.id, str(rec.seq)) for rec in chunk]
            futures.append(executor.submit(worker_process_chunk, lightweight, min_repeats, max_repeats, motive_size))

            for future in as_completed(futures):
                rows = future.result()
                if rows:
                    cur.executemany("INSERT OR IGNORE INTO repeats VALUES (?,?,?,?,?)", rows)
                    conn.commit()

    print("Writing results to output.txt...")
    write_repeats_to_txt(conn, "output.txt")
    conn.close()
    tmp.close()