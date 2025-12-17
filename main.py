from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from linkedList import Node, LinkedList

#fasta stückweise einlesen
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

#parallelisiertes einlesen
def process_sequence(seq):
    base_binary = {"A": "00", "C": "01", "G": "10", "T": "11"}
    base_tuple = make_dictionary()

    for i in range(len(seq) - 1):
        key = int(base_binary[seq[i]]+base_binary[seq[i + 1]],2)
        if base_tuple[key] is None:
            ll = LinkedList()
            base_tuple[key] = ll
            ll.append(i)
        else:
            base_tuple[key].append(i)

    return base_tuple

#test für hashmap(manuell)
def print_hashmap(base_tuple):
    index=0
    for _ in base_tuple:
        try:
            count = base_tuple[index].lenght()
        except:
            count = 0
        binary = bin(index)
        print(str(binary)+": "+str(count))
        index = index + 1


def make_dictionary():
    d = {}
    for i in range(0,16):
        d[i] = None
    return d


if __name__ == "__main__":
    fasta_path = "test_data_split.fasta"
    chunk_size = 5
    motive_size = 4
    min_repeats = 5
    max_repeats = 10

    #test
    chunk = fast_single_record(fasta_path)
    statistical_repeats(process_sequence(chunk[0].seq))

    #main

    # for batch in fasta_in_chunks(fasta_path, chunk_size):
    #
    #     with ProcessPoolExecutor(max_workers=len(batch)) as executor:
    #         results = list(executor.map(process_sequence, batch))
    #
    #     print("Batch fertig. Ergebnisse:")
    #     for r in results:
    #         print(r)