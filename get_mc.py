import multiprocessing
import argparse
import pysam
import collections
import subprocess
from Bio import SeqIO

def reverse_complement(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    reverse = dna[::-1].translate(complement)
    return reverse
def readgenome(name):
    dd = {}
    with open(name, 'r') as input_fasta:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            dd[record.id] = str(record.seq)
    return dd
def get_cut1_seq(fastq_file):
    seq_dict = {}

    with open(fastq_file) as f:
        while True:
            line = f.readline().strip()
            if not line:
                break  # 到达文件末尾，跳出循环
            seq_id = line.split()[0][1:]  # 去掉'@'字符
            seq = f.readline().strip()
            f.readline()  # 跳过不需要的行
            f.readline()  # 跳过不需要的行
            seq_dict[seq_id] = seq
    return seq_dict







def split_sam_file(sam_file, num_parts):
    with open(sam_file, "r") as f:
        lines = f.readlines()
    header_lines = []
    for line in lines:
        if line.startswith("@"):
            header_lines.append(line)
        else:
            break
    non_header_lines = lines[len(header_lines):]
    num_lines = len(non_header_lines)
    part_size = num_lines // num_parts
    parts = []
    for i in range(num_parts):
        start = i * part_size
        end = (i + 1) * part_size
        if i == num_parts - 1:
            end = None
        part_lines = header_lines + non_header_lines[start:end]
        parts.append(part_lines)
    for i, part in enumerate(parts):
        filename = f"{sam_file}.part{i + 1}.sam"
        with open(filename, "w") as f:
            for line in part:
                f.write(line)
    return [f'{sam_file}.part{i}.sam' for i in range(1, num_parts+1)], [f'{sam_file}.part{i}.sort' for i in range(1, num_parts+1)]



def process_file(bam_file):
    sam = pysam.AlignmentFile(bam_file, "r")
    hmc_bed = open(bam_file + ".bed", "w")
    for read in sam.fetch():
        if read.is_unmapped:
            continue
        cut_seq = ss[read.query_name]
        seq = read.seq
        aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=False)
        if read.is_reverse:
            cut_seq=reverse_complement(cut_seq)
            for query_pos, ref_pos, ref_base in aligned_pairs:
                if ref_pos is None or query_pos is None:
                    continue
                elif seq[query_pos] == "G" and cut_seq[query_pos] == "G" and dd[read.reference_name][ref_pos] == "G":
                    get_pos = read.reference_name + " " + str(ref_pos + 1) + " " + str(ref_pos + 1) + " "
                    if ref_pos == 0:
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif ref_pos == 1 and dd[read.reference_name][ref_pos - 1] != "C":
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 1] == "C":
                        print(get_pos + str(1) + " " + "CpG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 1] == "N" :
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 2] == "C":
                        print(get_pos + str(1) + " " + "CHG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 2] == "N" :
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    else:
                        print(get_pos + str(1) + " " + "CHH", file=hmc_bed)
                elif seq[query_pos] == "G" and cut_seq[query_pos] == "A" and ref_base == "G":
                    get_pos = read.reference_name + " " + str(ref_pos + 1) + " " + str(ref_pos + 1) + " "
                    if ref_pos == 0:
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif ref_pos == 1 and dd[read.reference_name][ref_pos - 1] != "C":
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 1] == "C":
                        print(get_pos + str(0) + " " + "CpG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 1] == "N" :
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 2] == "C":
                        print(get_pos + str(0) + " " + "CHG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 2] == "N" :
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    else:
                        print(get_pos + str(0) + " " + "CHH", file=hmc_bed)
        else:
            for query_pos, ref_pos, ref_base in aligned_pairs:
                if ref_pos is None or query_pos is None:
                    continue
                elif seq[query_pos] == "C" and cut_seq[query_pos] == "C" and dd[read.reference_name][ref_pos] == "C":
                    get_pos = read.reference_name + " " + str(ref_pos + 1) + " " + str(ref_pos + 1) + " "
                    if ref_pos == len(dd[read.reference_name])-1:
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif ref_pos == len(dd[read.reference_name])-2 and dd[read.reference_name][ref_pos+1]!="G":
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +1] == "G":
                        print(get_pos + str(1) + " " + "CpG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +1] == "N" :
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +2] == "G":
                        print(get_pos + str(1) + " " + "CHG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +2] == "N" :
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    else:
                        print(get_pos + str(1) + " " + "CHH", file=hmc_bed)
                elif seq[query_pos] == "C" and cut_seq[query_pos] == "T" and ref_base == "C":
                    get_pos = read.reference_name + " " + str(ref_pos + 1) + " " + str(ref_pos + 1) + " "
                    if ref_pos == len(dd[read.reference_name])-1:
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif ref_pos == len(dd[read.reference_name])-2 and dd[read.reference_name][ref_pos+1]!="G":
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +1] == "G":
                        print(get_pos + str(0) + " " + "CpG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +1] == "N" :
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +2] == "G":
                        print(get_pos + str(0) + " " + "CHG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +2] == "N" :
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    else:
                        print(get_pos + str(0) + " " + "CHH", file=hmc_bed)
    sam.close()
    hmc_bed.close()








def index_file(input_file):
    output_file = f"{input_file[:-4]}.sort"
    sort_command = f"samtools sort {input_file} > {output_file}"
    index_command = f"samtools index {output_file}"
    subprocess.run(sort_command, shell=True, check=True)
    subprocess.run(index_command, shell=True, check=True)



def final_bed(hmc_bed):
    my_bed = collections.defaultdict(lambda: (0, 0))
    with open(hmc_bed) as fin, open(hmc_bed+".merge", "w") as fout:
        for line in fin:
            line = line.strip().split()
            chrom, start, end, hmc, hmc_type = line[0], int(line[1]), int(line[2]), int(line[3]), line[4]
            key = (chrom, start, hmc_type)
            if hmc == 1:
                my_bed[key] = (my_bed[key][0] + 1, my_bed[key][1])
            elif hmc == 0:
                my_bed[key] = (my_bed[key][0], my_bed[key][1] + 1)

        # 将处理结果写入输出文件
        for key, counts in my_bed.items():
            chrom, start, hmc_type = key
            meth_ratio = counts[0] / (counts[0] + counts[1]) * 100
            meth_counts = counts[0]
            unmeth_counts = counts[1]
            fout.write(f"{chrom}\t{start}\t{start}\t{meth_ratio:.2f}\t{meth_counts}\t{unmeth_counts}\t{hmc_type}\n")


def process_chunk_chr(chunk):
    output_filename = args.sam+"." + chunk + ".bed"
    cmd = f"awk '$1==\"{chunk}\"' {output_file} > {output_filename}"
    subprocess.run(cmd, shell=True)




def run_command(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        print(stderr.decode())
    return stdout.decode()


def count_ccc(seq_name):
        seq = dd[seq_name]
        CpG_num = 0
        CHG_num = 0
        CHH_num = 0
        CN_num = 0
        for ref_pos in range(len(seq)):
            if seq[ref_pos] == 'C':
                if ref_pos == len(seq) - 1:
                    CN_num = CN_num + 1
                elif ref_pos == len(seq) - 2 and seq[ref_pos + 1] != 'G':
                    CN_num = CN_num + 1
                elif seq[ref_pos + 1] == 'G':
                    CpG_num = CpG_num + 1
                elif seq[ref_pos + 1] == 'N':
                    CN_num = CN_num + 1
                elif seq[ref_pos + 2] == 'G':
                    CHG_num = CHG_num + 1
                elif seq[ref_pos + 2] == 'N':
                    CN_num = CN_num + 1
                else:
                    CHH_num = CHH_num + 1
        return [CpG_num, CHG_num, CHH_num, CN_num]


def count_ggg(seq_name):
    seq = dd[seq_name]
    CpG_num = 0
    CHG_num = 0
    CHH_num = 0
    CN_num = 0
    for ref_pos in range(len(seq)):
        if seq[ref_pos] == 'G':
            if ref_pos == 0:
                CN_num = CN_num + 1
            elif ref_pos == 1 and seq[ref_pos + 1] != 'C':
                CN_num = CN_num + 1
            elif seq[ref_pos -1] == 'C':
                CpG_num = CpG_num + 1
            elif seq[ref_pos - 1] == 'N':
                CN_num = CN_num + 1
            elif seq[ref_pos - 2] == 'C':
                CHG_num = CHG_num + 1
            elif seq[ref_pos - 2] == 'N':
                CN_num = CN_num + 1
            else:
                CHH_num = CHH_num + 1
    return [CpG_num, CHG_num, CHH_num, CN_num]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--parallel",
        help="Number of parallel threads to use for processing. Default: 20",
        default=20,
        type=int,
    )

    parser.add_argument(
        "--ref",
        help="ref fasta",
        default="/data/wangzc/cgs/mc/new_mm9_lam_mc.fa",
        type=str,
    )

    parser.add_argument(
        "--sam",
        help="sam file",
        type=str,
    )

    parser.add_argument(
        "--output",
        help="output file",
        type=str,
    )

    parser.add_argument(
        "--cutfq1",
        help="output file",
        type=str,
    )

    args = parser.parse_args()
    if args.output is None:
        args.output = args.sam+".all.bed.mc"





    #读取参考基因组
    dd = readgenome(args.ref)
    ss = get_cut1_seq(args.cutfq1)
    input_split_files, sort_file = split_sam_file(args.sam, args.parallel)

    with multiprocessing.Pool(processes=args.parallel) as pool:
        pool.map(index_file, input_split_files)

    with multiprocessing.Pool(processes=args.parallel) as pool:
        pool.map(process_file, sort_file)


    #逐条提取已经结束，剩下全是整理合并






    output_file = args.sam+".merged.bed"
    command = f"cat {args.sam}.part*.sort.bed > {output_file}"
    subprocess.run(command, shell=True)

    delete_command = f"rm -rf {args.sam}.part*"
    subprocess.run(delete_command, shell=True)

    chunks_chr = [key for key in dd.keys()]

    with multiprocessing.Pool(processes=args.parallel) as pool:
        pool.map(process_chunk_chr, chunks_chr)

    delete_command = f"rm -rf {output_file}"
    subprocess.run(delete_command, shell=True)



    final_list=[args.sam+"." + key + ".bed" for key in chunks_chr]
    with multiprocessing.Pool(processes=args.parallel) as pool:
        pool.map(final_bed, final_list)

    delete_command = f"rm -rf {args.sam}.*.bed"
    subprocess.run(delete_command, shell=True)

    command = f"cat {args.sam}.*.bed.merge > {args.output}"
    subprocess.run(command, shell=True)

    delete_command = f"rm -rf {args.sam}.*.bed.merge"
    subprocess.run(delete_command, shell=True)






    with multiprocessing.Pool(processes=args.parallel) as pool:

        commands = [f"awk '$7==\"CpG\"{{sum4+=$4; sum5+=$5; sum6+=$6; cnt_cpg++}} END {{print cnt_cpg==0 ? 0 : sum4/cnt_cpg, sum5, sum6, cnt_cpg}}' {args.output}",
                f"awk '$7==\"CHG\"{{sum4+=$4; sum5+=$5; sum6+=$6; cnt_chg++}} END {{print cnt_chg==0 ? 0 : sum4/cnt_chg, sum5, sum6, cnt_chg}}' {args.output}",
                f"awk '$7==\"CHH\"{{sum4+=$4; sum5+=$5; sum6+=$6; cnt_chh++}} END {{print cnt_chh==0 ? 0 : sum4/cnt_chh, sum5, sum6, cnt_chh}}' {args.output}",
                f"awk '$7==\"CN\"{{sum4+=$4; sum5+=$5; sum6+=$6; cnt_cn++}} END {{print cnt_cn==0 ? 0 : sum4/cnt_cn, sum5, sum6, cnt_cn}}' {args.output}",
        ]

        ccc = pool.map(count_ccc, chunks_chr)
        ggg = pool.map(count_ggg, chunks_chr)
        results = pool.map(run_command, commands)

    ccc_CpG=ccc_CHG=ccc_CHH=ccc_CN=0
    for i in ccc:
        ccc_CpG=ccc_CpG+i[0]
        ccc_CHG=ccc_CHG+i[1]
        ccc_CHH=ccc_CHH+i[2]
        ccc_CN=ccc_CN+i[3]
    for i in ggg:
        ccc_CpG=ccc_CpG+i[0]
        ccc_CHG=ccc_CHG+i[1]
        ccc_CHH=ccc_CHH+i[2]
        ccc_CN=ccc_CN+i[3]

    print(results)

    result_cpg = int(results[0].strip().split()[3]) / ccc_CpG * 100
    result_cpg = round(result_cpg, 2)

    result_chg = int(results[1].strip().split()[3]) / ccc_CHG * 100
    result_chg = round(result_chg, 2)

    result_chh = int(results[2].strip().split()[3]) / ccc_CHH * 100
    result_chh = round(result_chh, 2)

    with open(args.output+".log","w")as file:
        print(f"CpG mc平均水平：{results[0].split()[0]} ref总CpG位点：{ccc_CpG} 测序覆盖到的CpG位点数：{results[0].strip().split()[3]} ({result_cpg}) modC总数量：{results[0].split()[1]} unmodC总数量：{results[0].split()[2]}",file=file)
        print(f"CHG mc平均水平：{results[1].split()[0]} ref总CHG位点：{ccc_CHG} 测序覆盖到的CHG位点数：{results[1].strip().split()[3]} ({result_chg}) modC总数量：{results[1].split()[1]} unmodC总数量：{results[1].split()[2]}",file=file)
        print(f"CHH mc平均水平：{results[2].split()[0]} ref总CHH位点：{ccc_CHH} 测序覆盖到的CHH位点数：{results[2].strip().split()[3]} ({result_chh}) modC总数量：{results[2].split()[1]} unmodC总数量：{results[2].split()[2]}",file=file)
        if float(results[3].split()[0]) != 0:
            print(f"CN mc平均水平：{results[3].split()[0]}   modC总数量：{results[3].strip().split()[1]} unmodC总数量：{results[3].split()[2]}",file=file)


















