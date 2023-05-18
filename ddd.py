with open('restored_seq.fq', 'r') as fin, \
        open('output1.fastq', 'w') as fout1, \
        open('output2.fastq', 'w') as fout2:
    line_num = 0
    for line in fin:
        line_num += 1
        if line_num % 4 == 1:  # header line
            header = line.strip()
        elif line_num % 4 == 2:  # sequence line
            seq = line.strip()
            mid = len(seq) // 2
            seq1, seq2 = seq[:mid], seq[mid:]
        elif line_num % 4 == 0:  # quality line
            qual = line.strip()
            fout1.write(f"{header}_1\n{seq1}\n+\n{qual[:mid]}\n")
            fout2.write(f"{header}_2\n{seq2}\n+\n{qual[mid:]}\n")

