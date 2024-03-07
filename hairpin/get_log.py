import argparse
from Bio import SeqIO
import numpy as np







def readgenome(name):
    dd = {}
    with open(name, 'r') as input_fasta:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            dd[record.id] = str(record.seq)
    return dd
def count_ref_C(dd):
    result=[]
    for seq in dd.values():
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

            elif seq[ref_pos] == 'G':
                if ref_pos == 0:
                    CN_num = CN_num + 1
                elif ref_pos == 1 and seq[ref_pos + 1] != 'C':
                    CN_num = CN_num + 1
                elif seq[ref_pos - 1] == 'C':
                    CpG_num = CpG_num + 1
                elif seq[ref_pos - 1] == 'N':
                    CN_num = CN_num + 1
                elif seq[ref_pos - 2] == 'C':
                    CHG_num = CHG_num + 1
                elif seq[ref_pos - 2] == 'N':
                    CN_num = CN_num + 1
                else:
                    CHH_num = CHH_num + 1

        my_string = ' '.join(map(str, [CpG_num, CHG_num, CHH_num, CN_num]))
        result.append(my_string)

    ccc_CpG = ccc_CHG = ccc_CHH = ccc_CN = 0
    for line in result:
        line = line.strip()
        line = line.split()
        ccc_CpG = ccc_CpG + int(line[0])
        ccc_CHG = ccc_CHG + int(line[1])
        ccc_CHH = ccc_CHH + int(line[2])
        ccc_CN = ccc_CN + int(line[3])


    return ccc_CpG,ccc_CHG,ccc_CHH,ccc_CN




def calculate_statistics(file_path,aaa,bbb,ccc,mode,spikein):
    cpg_sum4 = cpg_sum5 = cpg_sum6 = cpg_cnt = 0
    chg_sum4 = chg_sum5 = chg_sum6 = chg_cnt = 0
    chh_sum4 = chh_sum5 = chh_sum6 = chh_cnt = 0
    cn_sum4 = cn_sum5 = cn_sum6 = cn_cnt = 0


    sp_dir={}
    for i in spikein:
        sp_dir[i]=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    not_equal_list=set(spikein)
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split()
            if columns[0] not in not_equal_list:
                if columns[6] == 'CpG':
                    cpg_sum4 += float(columns[3])
                    cpg_sum5 += int(columns[4])
                    cpg_sum6 += int(columns[5])
                    cpg_cnt += 1
                elif columns[6] == 'CHG':
                    chg_sum4 += float(columns[3])
                    chg_sum5 += int(columns[4])
                    chg_sum6 += int(columns[5])
                    chg_cnt += 1
                elif columns[6] == 'CHH':
                    chh_sum4 += float(columns[3])
                    chh_sum5 += int(columns[4])
                    chh_sum6 += int(columns[5])
                    chh_cnt += 1
                else:
                    cn_sum4 += float(columns[3])
                    cn_sum5 += int(columns[4])
                    cn_sum6 += int(columns[5])
                    cn_cnt += 1
            else:
                if columns[6] == 'CpG':
                    sp_dir[columns[0]]=sp_dir[columns[0]] + np.array([float(columns[3]),1, 0.0, 0.0, 0.0, 0.0])
                if columns[6] == 'CHG':
                    sp_dir[columns[0]]=sp_dir[columns[0]] + np.array([0.0,0.0,float(columns[3]),1, 0.0, 0.0])
                if columns[6] == 'CHH':
                    sp_dir[columns[0]]=sp_dir[columns[0]] + np.array([0.0,0.0,0.0,0.0,float(columns[3]),1])


    avg_cpg = 100*cpg_sum4 / cpg_cnt if cpg_cnt != 0 else 0
    avg_chg = 100*chg_sum4 / chg_cnt if chg_cnt != 0 else 0
    avg_chh = 100*chh_sum4 / chh_cnt if chh_cnt != 0 else 0
    avg_cn =  100*cn_sum4 / cn_cnt if cn_cnt != 0 else 0

    if mode == "hmc":
        cpg = f"CpG hmc平均水平：{round(avg_cpg,2)}% ref总CpG位点：{aaa} 测序覆盖到的CpG位点数：{cpg_cnt} ({round(cpg_cnt*100/aaa, 2)}%) modC总数量：{cpg_sum5} unmodC总数量：{cpg_sum6}"
        chg = f"CHG hmc平均水平：{round(avg_chg,2)}% ref总CHG位点：{bbb} 测序覆盖到的CHG位点数：{chg_cnt} ({round(chg_cnt * 100 / bbb, 2)}%) modC总数量：{chg_sum5} unmodC总数量：{chg_sum6}"
        chh = f"CHH hmc平均水平：{round(avg_chh,2)}% ref总CHH位点：{ccc} 测序覆盖到的CHH位点数：{chh_cnt} ({round(chh_cnt * 100 / ccc, 2)}%) modC总数量：{chh_sum5} unmodC总数量：{chh_sum6}"
        cn = f"CN hmc平均水平：{round(avg_cn,2)}% modC总数量：{cn_sum5} unmodC总数量：{cn_sum6}"

        result=[cpg,chg,chh,cn]
        for key,value in sp_dir.items():
            avg_tempcpg = 100*value[0] / value[1] if value[1] != 0 else 0
            avg_tempchg = 100*value[2] / value[3] if value[3] != 0 else 0
            avg_tempchh = 100*value[4] / value[5] if value[5] != 0 else 0

            temp = f"hmc位点，{key}标品CpG位点假阳性率：{round(avg_tempcpg, 2)}% CHG位点假阳性率：{round(avg_tempchg, 2)}% CHH位点假阳性率：{round(avg_tempchh, 2)}%"
            result.append(temp)


    else :
        cpg = f"CpG mc平均水平：{round(avg_cpg,2)}% ref总CpG位点：{aaa} 测序覆盖到的CpG位点数：{cpg_cnt} ({round(cpg_cnt * 100 / aaa, 2)}%) modC总数量：{cpg_sum5} unmodC总数量：{cpg_sum6}"
        chg = f"CHG mc平均水平：{round(avg_chg,2)}% ref总CHG位点：{bbb} 测序覆盖到的CHG位点数：{chg_cnt} ({round(chg_cnt * 100 / bbb, 2)}%) modC总数量：{chg_sum5} unmodC总数量：{chg_sum6}"
        chh = f"CHH mc平均水平：{round(avg_chh,2)}% ref总CHH位点：{ccc} 测序覆盖到的CHH位点数：{chh_cnt} ({round(chh_cnt * 100 / ccc, 2)}%) modC总数量：{chh_sum5} unmodC总数量：{chh_sum6}"
        cn = f"CN mc平均水平：{round(avg_cn,2)}% modC总数量：{cn_sum5} unmodC总数量：{cn_sum6}"
        result = [cpg, chg, chh, cn]
        for key,value in sp_dir.items():
            avg_tempcpg = 100*value[0] / value[1] if value[1] != 0 else 0
            avg_tempchg = 100*value[2] / value[3] if value[3] != 0 else 0
            avg_tempchh = 100*value[4] / value[5] if value[5] != 0 else 0

            temp = f"mc位点，{key}标品CpG位点水平：{round(avg_tempcpg,2)}% CHG位点水平：{round(avg_tempchg,2)}% CHH位点水平：{round(avg_tempchh,2)}%"
            result.append(temp)


    return result




def main():
    parser = argparse.ArgumentParser(description="Duplicate Base Kmer Tool")

    parser.add_argument("--all_hmc", type=str, help="Input fq1 file")
    parser.add_argument("--all_mc", type=str, help="Input fq2 file")
    parser.add_argument("--spikein", type=str, help="default: 2") #逗号连接
    parser.add_argument("--output", type=str, help="default: 2")
    parser.add_argument("--ref", type=str, default="",help="Output w file")

    args = parser.parse_args()

    if args.ref=="":
        a = 43452244
        b = 229989601
        c = 819319197
    else:
        dd=readgenome(args.ref)
        a,b,c,d=count_ref_C(dd)


    spikein=args.spikein.split(',')
    hmc_list=calculate_statistics(args.all_hmc,a,b,c,"hmc",spikein)
    mc_list = calculate_statistics(args.all_mc, a, b, c, "mc", spikein)

    with open(args.output,"a")as file:
        for i in hmc_list:
            print(i,file=file)
        for i in mc_list:
            print(i,file=file)


if __name__ == '__main__':
    main()











