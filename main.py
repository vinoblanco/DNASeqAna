import itertools
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from linkedList import LinkedList

#fasta stückweise einlesen --> todo RAM gesteuert machen
def fasta_in_chunks(fasta_path, chunk_size):
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

#Hashmap erstellen und mit den linked lists befüllen
def process_sequence(seq):
    base_binary = {"A": "00", "C": "01", "G": "10", "T": "11"}
    #base_tuple = make_dictionary()
    base_tuple = {''.join(kombi): None for kombi in itertools.product(['A','C','G','T'], repeat=2)}
    print(base_tuple)

    for i in range(len(seq) - 1):
        #key = int(base_binary[seq[i]]+base_binary[seq[i + 1]],2)
        key = str(seq[i:i+2])
        if base_tuple[key] is None:
            ll = LinkedList()
            base_tuple[key] = ll
            ll.append(i)
        else:
            base_tuple[key].append(i)
    print_hashmap(base_tuple)

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

#liest die Hashmap ein und gibt die Repeats mit einigen Daten aus
def statistical_repeats(base_tuple, seq, min_repeats, max_repeats):
    repeats = []
    diffs = []
    print(base_tuple)

    for key, value in base_tuple.items():
        print(key)
        if value is None or value.lenght() < min_repeats:
            continue

        print(key)
        for counter, firststart in calculate_difference(value):
            print(counter)
            for period, count in counter.items():
                if count >= min_repeats - 1:
                    start = firststart
                    print(start)
                    end = start + period * count

                    if end <= len(seq):
                        motif = seq[start:start+period]

                        repeats.append({
                            "pair": key,
                            "start": start,
                            "period": period,
                            "repeat": count + 1,
                            "motif": motif
                        })

    return repeats


def make_dictionary():
    d = {}
    for i in range(0,16):
        d[i] = None
    return d

#Nimmt immer zwei Vorkommen von einem BasenTuple und vergleicht die Abstände, sobald sich dieser ändert wird es ausgegeben
def calculate_difference(ll):
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

if __name__ == "__main__":
    fasta_path = "test_data_split.fasta"
    chunk_size = 5
    motive_size = 4
    min_repeats = 5
    max_repeats = 10

    #test
    chunk = fast_single_record(fasta_path)
    print(statistical_repeats(process_sequence(chunk[0].seq), chunk[0], min_repeats, max_repeats))

    #main

    # for batch in fasta_in_chunks(fasta_path, chunk_size):
    #
    #     with ProcessPoolExecutor(max_workers=len(batch)) as executor:
    #         results = list(executor.map(process_sequence, batch))
    #
    #     print("Batch fertig. Ergebnisse:")
    #     for r in results:
    #         print(r)