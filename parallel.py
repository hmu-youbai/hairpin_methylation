import gzip
import math
from collections import namedtuple
import multiprocessing
from cgs.trim import trim_overlap
from functools import partial

Record = namedtuple('Record', ['name', 'seq', 'qual'])


def cut(read, qual):
    half_len = len(read)
    if half_len > len(qual):
        qual += 'I' * (half_len - len(qual))
    return qual[:half_len]


def process_chunk(chunks, alignment, min_lenth=50):
    result = []
    for record in chunks:
        read1, read2 = trim_overlap(record.seq.split()[0], record.seq.split()[1])
        seq = alignment.align(read1, read2)
        if len(seq) < min_lenth:
            continue
        qual = cut(seq, record.qual)
        result.append(Record(name=record.name, seq=seq, qual=qual))
    return result


def process_records(buffer, alignment, min_lenth=50, processes=24):
    pool = multiprocessing.Pool(processes=processes)
    chunk_size = math.ceil(len(buffer) / processes)
    chunks = [buffer[i:i + chunk_size] for i in range(0, len(buffer), chunk_size)]
    partial_process_chunk = partial(process_chunk, alignment=alignment, min_lenth=min_lenth)
    results = pool.map(partial_process_chunk, chunks)
    pool.close()
    pool.join()
    return [record for result in results for record in result]


def cache_and_process(input_file1, input_file2, output_file, alignment, min_lenth=50, processes=24, chunk_size=1000000):
    with gzip.open(input_file1, "rt") as in_handle, gzip.open(input_file2, "rt") as in_handle2, open(output_file,
                                                                                                     "w") as out_handle:
        buffer = []
        for line_num, (line1, line2) in enumerate(zip(in_handle, in_handle2)):
            if line_num % 4 == 0:
                record_name1 = line1.strip().split()[0]
            elif line_num % 4 == 1:
                record_seq = line1.strip() + " " + line2.strip()
            elif line_num % 4 == 3:
                record_qual1 = line1.strip()
                buffer.append(Record(name=record_name1, seq=record_seq, qual=record_qual1))
            if len(buffer) == chunk_size:
                print("process reads : " + str((line_num + 1) / 4))
                temp = process_records(buffer, alignment, min_lenth, processes)
                for seq in temp:
                    out_handle.write(f"{seq.name}\n")
                    out_handle.write(f"{seq.seq}\n")
                    out_handle.write(f"+\n")
                    out_handle.write(f"{seq.qual}\n")
                buffer = []
        if buffer:
            temp = process_records(buffer, alignment, min_lenth, processes)
            for seq in temp:
                out_handle.write(f"{seq.name}\n")
                out_handle.write(f"{seq.seq}\n")
                out_handle.write(f"+\n")
                out_handle.write(f"{seq.qual}\n")


