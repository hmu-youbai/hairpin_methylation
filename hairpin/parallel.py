import gzip
import math
from collections import namedtuple
import multiprocessing
from hairpin.trim import trim_overlap
from functools import partial
import os

Record = namedtuple('Record', ['name', 'seq', 'qual', 'seq1', 'seq2'])


def cut(read, qual):
    half_len = len(read)
    if half_len > len(qual):
        qual += 'I' * (half_len - len(qual))
    return qual[:half_len]


def process_chunk(chunks, alignment, min_lenth=50):
    result = []
    for record in chunks:
        read1, read2 = trim_overlap(record.seq.split()[0], record.seq.split()[1])
        seq, seq1, seq2 = alignment.align(read1, read2)
        if len(seq) < min_lenth:
            continue
        qual = cut(seq, record.qual)
        result.append(Record(name=record.name, seq=seq, qual=qual, seq1=seq1, seq2=seq2))
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
    # 获取output_file所在的目录
    output_dir = os.path.dirname(output_file)

    # 创建cut_f1.fq和cut_f2.fq文件
    cut_f1_path = os.path.join(output_dir, output_file[:-2]+'cut_f1.fq')
    cut_f2_path = os.path.join(output_dir, output_file[:-2]+'cut_f2.fq')
    report_n1 = 0
    if input_file1.endswith('.gz'):
        in_handle = gzip.open(input_file1, "rt")
    else:
        in_handle = open(input_file1)
    if input_file2.endswith('.gz'):
        in_handle2 = gzip.open(input_file2, "rt")
    else:
        in_handle2 = open(input_file2)
    with in_handle, in_handle2, open(output_file, "w") as out_handle, open(cut_f1_path, "w") as f1, open(cut_f2_path, "w") as f2:
        buffer = []
        for line_num, (line1, line2) in enumerate(zip(in_handle, in_handle2)):
            if line_num % 4 == 0:
                record_name1 = line1.strip().split()[0]
            elif line_num % 4 == 1:
                record_seq = line1.strip() + " " + line2.strip()
            elif line_num % 4 == 3:
                record_qual1 = line1.strip()
                buffer.append(Record(name=record_name1, seq=record_seq, qual=record_qual1, seq1=line1.strip(),
                                     seq2=line2.strip()))
            if len(buffer) == chunk_size:
                print("process reads : " + str(int((line_num + 1) / 4)))
                temp = process_records(buffer, alignment, min_lenth, processes)
                report_n1 = report_n1 + len(temp)
                for seq in temp:
                    out_handle.write(f"{seq.name}\n")
                    out_handle.write(f"{seq.seq}\n")
                    out_handle.write(f"+\n")
                    out_handle.write(f"{seq.qual}\n")
                    f1.write(f"{seq.name}\n")
                    f1.write(f"{seq.seq1}\n")
                    f1.write(f"+\n")
                    f1.write(f"{seq.qual}\n")
                    f2.write(f"{seq.name}\n")
                    f2.write(f"{seq.seq2}\n")
                    f2.write(f"+\n")
                    f2.write(f"{seq.qual}\n")
                buffer = []
        if buffer:
            temp = process_records(buffer, alignment, min_lenth, processes)
            report_n1 = report_n1 + len(temp)
            for seq in temp:
                out_handle.write(f"{seq.name}\n")
                out_handle.write(f"{seq.seq}\n")
                out_handle.write(f"+\n")
                out_handle.write(f"{seq.qual}\n")
                f1.write(f"{seq.name}\n")
                f1.write(f"{seq.seq1}\n")
                f1.write(f"+\n")
                f1.write(f"{seq.qual}\n")
                f2.write(f"{seq.name}\n")
                f2.write(f"{seq.seq2}\n")
                f2.write(f"+\n")
                f2.write(f"{seq.qual}\n")
    in_handle.close()
    in_handle2.close()
    message = "fq1: " + input_file1 + "\n" + "fq2: " + input_file2 + "\n" + "min_read_length: " + str(
        min_lenth) + "\n" + "input_reads: " + str(int((line_num + 1) / 4)) + "\n" + "resolved_reads: " + str(
        report_n1) + "\n" + "resolved ratio: " + str("{:.2f}%".format(report_n1 / int((line_num + 1) / 4) * 100))
    print(message)
    return message
