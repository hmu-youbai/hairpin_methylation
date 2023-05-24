import numpy as np
import pandas as pd
from Bio import SeqIO


def get_lengths(filename):
    # 初始化序列长度列表
    lengths = []

    # 遍历fastq文件中的所有序列，并将它们的长度添加到列表中
    for record in SeqIO.parse(filename, "fastq"):
        lengths.append(len(record.seq))

    # 计算序列长度的统计信息
    mean_length = np.mean(lengths)
    median_length = np.median(lengths)
    std_length = np.std(lengths)

    # 将序列长度列表转换为Pandas Series对象
    lengths_series = pd.Series(lengths)

    # 计算每个长度的数量和占比
    counts = lengths_series.value_counts(normalize=True)

    # 获取占比最多的序列长度及其所占比例
    most_common_length = counts.index[0]
    most_common_length_percentage = counts.iloc[0]

    result_str = "Mean length: " + str(mean_length) + "\n" + \
                 "Median length: " + str(median_length) + "\n" + \
                 "Standard deviation: " + str(std_length) + "\n" + \
                 "Most common length: " + str(most_common_length) + "\n" + \
                 "Most common length percentage: " + str(most_common_length_percentage) + "\n"

    return result_str


def cover_for_log(mylog):
    cover_all = []
    cover_log = []
    dep = []
    mode = 0
    mc_log = 0
    mc_dep = 0
    with open(mylog) as file:
        for line in file:
            line = line.strip()
            line = line.split()

            if len(line) == 0:
                cover_all = cover_log
                cover_log = []
                lambda_all = lambda_log
                mc_all = mc_log
                dep = []

                my_lambda = lambda_dep
                my_mc = mc_dep
            elif len(line) == 11:
                if line[1] != "N":
                    cover_log.append(int(line[1].replace(',', '')))
                    dep.append(float(line[2].replace(',', '')))
                if line[0] == "lambda":
                    lambda_log = int(line[1].replace(',', ''))
                    lambda_dep = float(line[2].replace(',', ''))
                if line[0] == "mc":
                    mc_log = int(line[1].replace(',', ''))
                    mc_dep = float(line[2].replace(',', ''))
                    mode = 1

    ss = [dep[i] * cover_log[i] for i in range(len(dep))]

    if mode == 1:
        return sum(cover_log) / sum(cover_all), lambda_log / lambda_all, mc_log / mc_all, sum(ss) / sum(
            cover_log), my_lambda, my_mc
    else:
        return sum(cover_log) / sum(cover_all), lambda_log / lambda_all, 0, sum(ss) / sum(cover_log), my_lambda, my_mc
# cover,lambda_cover,mc_cover,dep,my_lambda,my_mc=cover_for_log("c2t.log")


def calculate_statistics(file_path,aaa,bbb,ccc,mode,args):
    cpg_sum4 = cpg_sum5 = cpg_sum6 = cpg_cnt = 0
    chg_sum4 = chg_sum5 = chg_sum6 = chg_cnt = 0
    chh_sum4 = chh_sum5 = chh_sum6 = chh_cnt = 0
    cn_sum4 = cn_sum5 = cn_sum6 = cn_cnt = 0

    mc_cpg=mc_chg=mc_chh=mccpg_cnt=mcchg_cnt=mcchh_cnt=0
    lambda_cpg = lambda_chg = lambda_chh = lambdacpg_cnt = lambdachg_cnt = lambdachh_cnt = 0
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split()
            if columns[0]!="mc" and  columns[0]!="lambda":
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

            elif columns[0]=="mc":
                if columns[6] == 'CpG':
                    mc_cpg += float(columns[3])
                    mccpg_cnt += 1
                elif columns[6] == 'CHG':
                    mc_chg += float(columns[3])
                    mcchg_cnt += 1
                elif columns[6] == 'CHH':
                    mc_chh += float(columns[3])
                    mcchh_cnt += 1

            elif columns[0]=="lambda":
                if columns[6] == 'CpG':
                    lambda_cpg += float(columns[3])
                    lambdacpg_cnt += 1
                elif columns[6] == 'CHG':
                    lambda_chg += float(columns[3])
                    lambdachg_cnt += 1
                elif columns[6] == 'CHH':
                    lambda_chh += float(columns[3])
                    lambdachh_cnt += 1

    avg_cpg = cpg_sum4 / cpg_cnt if cpg_cnt != 0 else 0
    avg_chg = chg_sum4 / chg_cnt if chg_cnt != 0 else 0
    avg_chh = chh_sum4 / chh_cnt if chh_cnt != 0 else 0
    avg_cn = cn_sum4 / cn_cnt if cn_cnt != 0 else 0

    avg_mccpg = mc_cpg / mccpg_cnt if mccpg_cnt != 0 else 0
    avg_mcchg = mc_chg / mcchg_cnt if mcchg_cnt != 0 else 0
    avg_mcchh = mc_chh / mcchh_cnt if mcchh_cnt != 0 else 0

    avg_lambdacpg = lambda_cpg / lambdacpg_cnt if lambdacpg_cnt != 0 else 0
    avg_lambdachg = lambda_chg / lambdachg_cnt if lambdachg_cnt != 0 else 0
    avg_lambdachh = lambda_chh / lambdachh_cnt if lambdachh_cnt != 0 else 0
    if mode == "hmc":
        cpg = f"CpG hmc平均水平：{avg_cpg} ref总CpG位点：{aaa} 测序覆盖到的CpG位点数：{cpg_cnt} ({round(cpg_cnt*100/aaa, 2)}) modC总数量：{cpg_sum5} unmodC总数量：{cpg_sum6}"
        chg = f"CHG hmc平均水平：{avg_chg} ref总CpG位点：{bbb} 测序覆盖到的CpG位点数：{chg_cnt} ({round(chg_cnt * 100 / bbb, 2)}) modC总数量：{chg_sum5} unmodC总数量：{chg_sum6}"
        chh = f"CHH hmc平均水平：{avg_chh} ref总CpG位点：{ccc} 测序覆盖到的CpG位点数：{chh_cnt} ({round(chh_cnt * 100 / ccc, 2)}) modC总数量：{chh_sum5} unmodC总数量：{chh_sum6}"
        cn = f"CN hmc平均水平：{avg_cn} modC总数量：{cn_sum5} unmodC总数量：{cn_sum6}"

        mc= f"mc标品CpG位点假阳性率：{avg_mccpg} CHG位点假阳性率：{avg_mcchg} CHH位点假阳性率：{avg_mcchh}"
        la = f"lambda标品CpG位点假阳性率：{avg_lambdacpg} CHG位点假阳性率：{avg_lambdachg} CHH位点假阳性率：{avg_lambdachh}"
    else :
        cpg = f"CpG mc平均水平：{avg_cpg} ref总CpG位点：{aaa} 测序覆盖到的CpG位点数：{cpg_cnt} ({round(cpg_cnt * 100 / aaa, 2)}) modC总数量：{cpg_sum5} unmodC总数量：{cpg_sum6}"
        chg = f"CHG mc平均水平：{avg_chg} ref总CpG位点：{bbb} 测序覆盖到的CpG位点数：{chg_cnt} ({round(chg_cnt * 100 / bbb, 2)}) modC总数量：{chg_sum5} unmodC总数量：{chg_sum6}"
        chh = f"CHH mc平均水平：{avg_chh} ref总CpG位点：{ccc} 测序覆盖到的CpG位点数：{chh_cnt} ({round(chh_cnt * 100 / ccc, 2)}) modC总数量：{chh_sum5} unmodC总数量：{chh_sum6}"
        cn = f"CN mc平均水平：{avg_cn} modC总数量：{cn_sum5} unmodC总数量：{cn_sum6}"

        mc = f"mc标品脱甲基CpG位点假阴性率：{100-avg_mccpg} CHG位点假阴性率：{100-avg_mcchg} CHH位点假阴性率：{100-avg_mcchh}"
        la = f"lambda标品CpG位点转化率：{100-avg_lambdacpg} CHG位点转化率：{100-avg_lambdachg} CHH位点转化率：{100-avg_lambdachh}"

    with open(args.output + ".log", "w") as file:
        print(cpg,file=file)
        print(chg, file=file)
        print(chh, file=file)
        print(cn, file=file)
        print(mc, file=file)
        print(la, file=file)












