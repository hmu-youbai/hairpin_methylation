import argparse
from hairpin.align import Alignment
from hairpin.parallel import cache_and_process
from hairpin.trim import deduplication
from hairpin.log import get_lengths
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
    mylog=mylog+"\n"+get_lengths(args.output)
    if args.log == 1:
        with open(args.fq1.rsplit(".", 1)[0] + ".log", "w") as file:
            print(mylog, file=file)




if __name__ == '__main__':
    main()
