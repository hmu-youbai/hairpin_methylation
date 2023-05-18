import gzip
from collections import namedtuple




Record = namedtuple('Record', ['name', 'qual'])
de_seq={}
with gzip.open("CGS-1w-1th_1.clean.fq.gz", "rt") as handle1,gzip.open("CGS-1w-1th_2.clean.fq.gz", "rt") as handle2:
    for line_num, (line1, line2) in enumerate(zip(handle1, handle2)):
        if line_num % 4 == 0:
            record_name1 = line1.split()[0]
        elif line_num % 4 == 1:
            record_seq = line1.strip() + " " + line2.strip()
        elif line_num % 4 == 3:
            de_seq[record_seq] = Record(name=record_name1,  qual=line1.strip() + " " + line2.strip())
    print(1111)

    with open("uniq1.fq", "w") as outfile1,open("uniq2.fq.gz", "w") as outfile2:
        for seq in de_seq.keys():
            print(de_seq[seq].name,file=outfile1)
            print(seq.split()[0],file=outfile1)
            print("+",file=outfile1)
            print(de_seq[seq].qual.split()[0],file=outfile1)
            print(de_seq[seq].name,file=outfile2)
            print(seq.split()[1],file=outfile2)
            print("+",file=outfile2)
            print(de_seq[seq].qual.split()[1],file=outfile2)



