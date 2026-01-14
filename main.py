import itertools
import sqlite3, tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from linkedList import LinkedList

#todo motive abgleichen
#todo optimieren das nicht komplette Sequenz mehrfach analysiert

#fasta stückweise einlesen --> todo RAM gesteuert machen
def fasta_in_chunks(fasta_path: str, chunk_size):
    """
     test
    :param fasta_path: Pfad zur Fasta Datei
    :param chunk_size: Wie viele Sequenzen pro Chunk berechnet werden
    :return: generiert Chunks mit Sequenzen
    """
    chunk = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        chunk.append(record)
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
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

def process_sequence(seq):
    """
    In dieser Methode wird ein Dictionary erstellt, welches mit Basentupeln als Keys und den vorkommen dieser als
    Value befüllt wird.
    :param seq: die Sequenz, für welche das Dic erstellt wird
    :return: das Dictionary
    """
    base_tuple = {''.join(kombi): None for kombi in itertools.product(['A','C','G','T'], repeat=2)}
    #print(base_tuple)

    for i in range(len(seq) - 1):
        key = str(seq[i:i+2])
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
        print(key)
        try:
            count = value.lenght()
            value.iterate()
        except:
            count = 0

def statistical_repeats(base_tuple, seq, min_repeats, max_repeats):
    """
    Werdet das Dictionary aus und schreibt die gefundenen Repeats in eine Datenbank
    :param base_tuple: das Dictionary
    :param seq: die Sequenz des Dictionary
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

                    if end <= len(seq):
                        id = seq.id
                        motif = str(canonical_dna_motif(seq[start:start+period].seq))
                        #print(motif)
                        repeats.append((id, motif, start, end, count))

    cur.executemany(
        "INSERT OR IGNORE INTO repeats VALUES (?,?,?,?,?)",
        repeats
    )
    conn.commit()

#todo motive vergleichen (hash vergleichen)
def calculate_difference(ll):
    """
    Wertet die Liked List aus, welche alle vorkommen eines Basentupels enthalten. Berechnet ob der abstand zwischen
    zwei Vorkommen periodisch ist.
    :param ll: Linked List mit den Vorkommen
    :return: Wann das Vorkommen startet (firststart) und wie oft es auftritt (count)
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
    Sortiert eine Sequenz lexiographisch und gibt das kleinste Ergebnis zurück
    :param seq: Sequenz
    :return: kleinste Sortierung
    """
    candidates = [seq[i:] + seq[:i] for i in range(len(seq))]
    return min(candidates)

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
        seq_number integer NOT NULL,
        motif text NOT NULL,
        start integer NOT NULL,
        period integer NOT NULL,
        repeat integer NOT NULL,
        UNIQUE (seq_number, motif))
    """)

    cur.execute("CREATE INDEX idx_motif ON repeats (motif)")
    conn.commit()

    #test
    chunk = fast_single_record(fasta_path)
    statistical_repeats(process_sequence(chunk[0].seq), chunk[0], min_repeats, max_repeats)
    cur.execute("""SELECT * FROM repeats ORDER BY motif""")
    for row in cur.fetchall():
        print(row)

    #main

    # for batch in fasta_in_chunks(fasta_path, chunk_size):
    #
    #     with ProcessPoolExecutor(max_workers=len(batch)) as executor:
    #         results = list(executor.map(process_sequence, batch))
    #
    #     print("Batch fertig. Ergebnisse:")
    #     for r in results:
    #         print(r)