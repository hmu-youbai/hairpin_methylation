import multiprocessing
import argparse
import pysam
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
    return [f'{sam_file}.part{i}.sam' for i in range(1, num_parts + 1)], [f'{sam_file}.part{i}.sort' for i in
                                                                          range(1, num_parts + 1)], header_lines


def process_file(bam_file):
    sam = pysam.AlignmentFile(bam_file, "r")
    sam_bis = open(bam_file[:-4] + "bis.sam", "w")
    for read in sam.fetch():
        if read.is_unmapped:
            continue
        cut_seq = ss[read.query_name]
        seq = read.seq
        aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=False)

        XM = []
        if read.is_reverse:
            XR = "XR:Z:CT\tXG:Z:GA"
            cut_seq = reverse_complement(cut_seq)
            for query_pos, ref_pos, ref_base in aligned_pairs:
                if ref_pos is None:
                    XM.append(".")
                elif query_pos is None:
                    continue
                elif seq[query_pos]=="G" and cut_seq[query_pos] == "G" and dd[read.reference_name][ref_pos] == "G":
                    if ref_pos == 0:
                        XM.append("U")
                    elif ref_pos == 1 and dd[read.reference_name][ref_pos - 1] != "C":
                        XM.append("U")
                    elif dd[read.reference_name][ref_pos - 1] == "C":
                        XM.append("Z")
                    elif dd[read.reference_name][ref_pos - 1] == "N":
                        XM.append("U")
                    elif dd[read.reference_name][ref_pos - 2] == "C":
                        XM.append("X")
                    elif dd[read.reference_name][ref_pos - 2] == "N":
                        XM.append("U")
                    else:
                        XM.append("H")
                elif seq[query_pos]=="G" and  cut_seq[query_pos] == "A" and ref_base == "G":
                    if ref_pos == 0:
                        XM.append("u")
                    elif ref_pos == 1 and dd[read.reference_name][ref_pos - 1] != "C":
                        XM.append("u")
                    elif dd[read.reference_name][ref_pos - 1] == "C":
                        XM.append("z")
                    elif dd[read.reference_name][ref_pos - 1] == "N":
                        XM.append("u")
                    elif dd[read.reference_name][ref_pos - 2] == "C":
                        XM.append("x")
                    elif dd[read.reference_name][ref_pos - 2] == "N":
                        XM.append("u")
                    else:
                        XM.append("h")
                else:
                    XM.append(".")

        else:
            XR = "XR:Z:CT\tXG:Z:CT"
            for query_pos, ref_pos, ref_base in aligned_pairs:
                if ref_pos is None:
                    XM.append(".")
                elif query_pos is None:
                    continue
                elif seq[query_pos]=="C" and cut_seq[query_pos] == "C" and dd[read.reference_name][ref_pos] == "C":
                    if ref_pos == len(dd[read.reference_name])-1:
                        XM.append("U")
                    elif ref_pos == len(dd[read.reference_name])-2 and dd[read.reference_name][ref_pos+1]!="G":
                        XM.append("U")
                    elif dd[read.reference_name][ref_pos +1] == "G":
                        XM.append("Z")
                    elif dd[read.reference_name][ref_pos +1] == "N" :
                        XM.append("U")
                    elif dd[read.reference_name][ref_pos +2] == "G":
                        XM.append("X")
                    elif dd[read.reference_name][ref_pos +2] == "N" :
                        XM.append("U")
                    else:
                        XM.append("H")
                elif seq[query_pos]=="C" and cut_seq[query_pos] == "T" and ref_base == "C":
                    if ref_pos == len(dd[read.reference_name])-1:
                        XM.append("u")
                    elif ref_pos == len(dd[read.reference_name])-2 and dd[read.reference_name][ref_pos+1]!="G":
                        XM.append("u")
                    elif dd[read.reference_name][ref_pos +1] == "G":
                        XM.append("z")
                    elif dd[read.reference_name][ref_pos +1] == "N" :
                        XM.append("u")
                    elif dd[read.reference_name][ref_pos +2] == "G":
                        XM.append("x")
                    elif dd[read.reference_name][ref_pos +2] == "N" :
                        XM.append("u")
                    else:
                        XM.append("h")
                else:
                    XM.append(".")

        nm = 0
        temp_num = 0
        past_query_pos = ""
        MD = ""
        for query_pos, ref_pos, ref_base in aligned_pairs:
            if query_pos is None:
                nm = nm + 1
                if past_query_pos is None:
                    MD = MD + dd[read.reference_name][ref_pos]
                else:
                    MD = MD + str(temp_num) + "^" + dd[read.reference_name][ref_pos]
                    temp_num = 0
            elif ref_pos is None:
                nm = nm + 1
            elif cut_seq[query_pos] == dd[read.reference_name][ref_pos]:
                temp_num = temp_num + 1
            else:
                nm = nm + 1
                MD = MD + str(temp_num) + dd[read.reference_name][ref_pos]
                temp_num = 0
            past_query_pos = query_pos
        MD = MD + str(temp_num)

        sam_result = read.query_name + "\t" + str(read.flag) + "\t" + read.reference_name + "\t" + str(
            read.reference_start+1) + "\t" + str(
            read.mapq) + "\t" + read.cigarstring + "\t" + "*" + "\t" + "0" + "\t" + "0" + "\t" + cut_seq + "\t" + read.qual + "\t" + "NM:i:" + str(
            nm) + "\t" + "MD:Z:" + MD + "\t" + "XM:Z:" + ''.join(XM) + "\t" + XR
        print(sam_result, file=sam_bis)
    sam.close()
    sam_bis.close()


def index_file(input_file):
    output_file = f"{input_file[:-4]}.sort"
    sort_command = f"samtools sort {input_file} > {output_file}"
    index_command = f"samtools index {output_file}"
    subprocess.run(sort_command, shell=True, check=True)
    subprocess.run(index_command, shell=True, check=True)




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
        args.output = args.sam + ".mc.all.bed"

    # 读取参考基因组
    dd = readgenome(args.ref)
    ss = get_cut1_seq(args.cutfq1)
    input_split_files, sort_file,header_lines = split_sam_file(args.sam, args.parallel)

    with multiprocessing.Pool(processes=args.parallel) as pool:
        pool.map(index_file, input_split_files)

    with multiprocessing.Pool(processes=args.parallel) as pool:
        pool.map(process_file, sort_file)


    output_file = args.sam+".merged.bis.sam"
    command = f"cat {args.sam}.part*.bis.sam > {output_file}"
    subprocess.run(command, shell=True)

    delete_command = f"rm -rf {args.sam}.part*"
    subprocess.run(delete_command, shell=True)



    with open(output_file, 'r') as f:
        original_content = f.read()
    new_content = ''.join(header_lines) + original_content
    with open(output_file, 'w') as f:
        f.write(new_content)





