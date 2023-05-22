import argparse
import pysam
from collections import Counter

def reverse_complement(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    reverse = dna[::-1].translate(complement)
    return reverse
def process_file(bam_file):
    sam = pysam.AlignmentFile(bam_file, "r")
    hmc_bed = open(bam_file + ".fq", "w")

    aaa=0

    for read in sam.fetch():
        aaa=aaa+1
        if aaa%10000==0:
            print(aaa)
        if read.is_unmapped:
            continue
        name=read.query_name

        aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=False)
        mylist = [pair[2] for pair in aligned_pairs if pair[0] is not None]
        new_list = [element.upper() if isinstance(element, str) else 'N' if element is None else element for element in mylist]

        mystr=''.join(new_list)
        if read.is_reverse:
            mystr=reverse_complement(mystr)

        print(name,file=hmc_bed)
        print(mystr,file=hmc_bed)
    sam.close()
    hmc_bed.close()







if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fq1",
        help="输入fq1"
    )
    parser.add_argument(
        "--fq2",
        help="输入fq2"
    )

    parser.add_argument(
        "--sam",
        default="restored_seq.fq.sam",
        help="输入sam"
    )



    args = parser.parse_args()
    process_file(args.sam)

    cc = {}
    with open(args.sam+".fq") as file:
        for i, line in enumerate(file):
            if i % 2 == 0:
                name = line.strip()
            else:
                cc[name] = line.strip()

    reff = []

    fq1_list = []
    with open(args.fq1, 'r') as f1:
        for i, line in enumerate(f1):
            if i % 4 == 0:
                name = line.strip()[1:]
            elif i % 4 == 1:
                if name in cc.keys():
                    fq1_list.append(line.strip())
                    reff.append(cc[name])

    fq2_list = []
    with open(args.fq2, 'r') as f2:
        for i, line in enumerate(f2):
            if i % 4 == 0:
                name = line.strip()[1:]
            elif i % 4 == 1:
                if name in cc.keys():
                    fq2_list.append(line.strip())

    print("fq1序列数目：", len(fq1_list))
    print("fq2序列数目：", len(fq2_list))

    tt = Counter(
        j + k + l for read1, read2, read3 in zip(fq1_list, fq2_list, reff) for j, k, l in zip(read1, read2, read3))

    counts = tt

    total_counts = sum(counts.values())
    with open(args.sam+".log","w")as file:
        for key, value in sorted(counts.items(), key=lambda x: x[1], reverse=True):
            percentage = '{:.2f}%'.format(value / total_counts * 100)
            print(key, value, percentage,file=file)























































