import argparse

from hairpin.align import Alignment
from hairpin.parallel import cache_and_process
from hairpin.trim import deduplication
from hairpin.log import get_lengths
from hairpin.log import cover_for_log
import subprocess
import os


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fq1",
        help="Input filename for forward read sequences.",
    )

    parser.add_argument(
        "--fq2",
        help="Input filename for reverse read sequences.",
    )

    parser.add_argument(
        "--output",
        help="Output filename for processed sequences. Default: restored_seq.fq",
    )

    parser.add_argument(
        "--rule",
        help="Choose a rule for sequence alignment: 1 or 2. Default: 1",
        default="2",
        choices=["1", "2"],
        type=str,
    )

    parser.add_argument(
        "--min",
        help="Minimum length of a sequence to be processed. Default: 50",
        default=50,
        type=int,
    )

    parser.add_argument(
        "--parallel",
        help="Number of parallel threads to use for processing. Default: 24",
        default=24,
        type=int,
    )

    parser.add_argument(
        "--chunk_size",
        help="Chunk size (number of sequences) for processing. Default: 10000000",
        default=10000000,
        type=int,
    )

    parser.add_argument(
        "--duplication",
        help="using deduplication",
        default=0,
        type=int,
    )

    parser.add_argument(
        "--rm_temp",
        help="rm temp file ",
        default=1,
        type=int,
    )

    parser.add_argument(
        "--log",
        help="get log",
        default=1,
        type=int,
    )

    parser.add_argument(
        "--bowtie2_ref",
        help="get log",
        default="/data/wangzc/cgs/mc/mm9",
        type=str,
    )

    parser.add_argument(
        "--BAMStats",
        help="get log",
        default="/data/wangzc/soft/BAMStats-1.25/BAMStats-1.25.jar",
        type=str,
    )

    args = parser.parse_args()

    if args.output is None:
        output_dir = os.path.dirname(args.fq1)
        args.output = os.path.join(output_dir, 'restored_seq.fq')

    if args.rule == "1":
        alignment = Alignment("rule_matrix_1")
    elif args.rule == "2":
        alignment = Alignment("rule_matrix_2")
    else:
        raise ValueError("Invalid output value. Valid values are '1' or '2'.")

    # 处理文件
    if args.duplication == 0:
        mylog = cache_and_process(args.fq1, args.fq2, args.output, alignment, args.min, args.parallel, args.chunk_size)
    else:
        uniq_read1, uniq_read2, raw_read = deduplication(args.fq1, args.fq2)
        print("read before deduplication: " + str(raw_read))
        mylog_1 = cache_and_process(uniq_read1, uniq_read2, args.output, alignment, args.min, args.parallel,
                                    args.chunk_size)
        mylog = "read before deduplication: " + str(raw_read) + "\n" + mylog_1
        if args.rm_temp == 1:
            cmd = ["rm", uniq_read1, uniq_read2]
            subprocess.call(cmd)
    mylog = mylog +"\n"+ get_lengths(args.output)


    bowtie2_cmd = f"bowtie2 -x  {args.bowtie2_ref} -U {args.output} -S {args.output}.sam -p {args.parallel}"
    mylog_bowtie2 = subprocess.run(bowtie2_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    mylog = mylog + "\n" + mylog_bowtie2.stdout.decode('utf-8')
    mylog = mylog + "\n" + mylog_bowtie2.stderr.decode('utf-8')


    bam_cmd = f"samtools view  -b {args.output}.sam > {args.output}.bam"
    sort_cmd= f"samtools sort {args.output}.bam > {args.output}.sort"
    index_cmd= f"samtools index  {args.output}.sort"
    cover_cmd= f"java -jar -Xmx100g {args.BAMStats} -m -i {args.output}.sort -o {args.output}.BAMStats --view simple"
    subprocess.run(bam_cmd, shell=True)
    subprocess.run(sort_cmd, shell=True)
    subprocess.run(index_cmd, shell=True)
    subprocess.run(cover_cmd, shell=True)
    cover, lambda_cover, mc_cover, dep, my_lambda, my_mc = cover_for_log(args.output+".BAMStats")
    mylog_len="基因组覆盖度为: " + str(cover) + "\n" + \
                 "lambda cover: " + str(lambda_cover) + "\n" + \
                 "full mc 标品 cover: " + str(mc_cover) + "\n" + \
                 "平均测序深度为: " + str(dep) + "\n" + \
                 "lambda测序深度为: " + str(my_lambda) + "\n" + \
                 "mc标品测序深度为: " + str(my_mc) + "\n"

    mylog = mylog +"\n" +mylog_len
    hmc_cmd=f"hmc_extractor  --sam {args.output}.sam --parallel {args.parallel}"
    mc_cmd=f"mc_extractor  --sam {args.output}.sam  --cutfq1 {args.output[:-2]}cut_f1.fq  --parallel {args.parallel}"
    subprocess.run(hmc_cmd, shell=True)
    subprocess.run(mc_cmd, shell=True)




    if args.log == 1:
        with open(args.fq1.rsplit(".", 1)[0] + ".log", "w") as file:
            print(mylog, file=file)





if __name__ == '__main__':
    main()
