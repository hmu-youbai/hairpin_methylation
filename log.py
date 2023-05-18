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
    dep=[]
    mode=0
    mc_log=0
    mc_dep=0
    with open(mylog) as file:
        for line in file:
            line = line.strip()
            line = line.split()

            if len(line) == 0:
                cover_all = cover_log
                cover_log = []
                lambda_all=lambda_log
                mc_all = mc_log
                dep=[]

                my_lambda=lambda_dep
                my_mc=mc_dep
            elif len(line) == 11:
                if line[1] != "N":
                    cover_log.append(int(line[1].replace(',', '')))
                    dep.append(float(line[2].replace(',', '')))
                if line[0]=="lambda":
                    lambda_log = int(line[1].replace(',', ''))
                    lambda_dep= float(line[2].replace(',', ''))
                if line[0]=="mc":
                    mc_log = int(line[1].replace(',', ''))
                    mc_dep= float(line[2].replace(',', ''))
                    mode=1

    ss = [dep[i] * cover_log[i] for i in range(len(dep))]

    if mode==1:
        return sum(cover_log)/sum(cover_all), lambda_log/lambda_all, mc_log/mc_all,sum(ss)/sum(cover_log),my_lambda,my_mc
    else:
        return sum(cover_log) / sum(cover_all), lambda_log / lambda_all, 0, sum(ss) / sum(cover_log), my_lambda, my_mc
# cover,lambda_cover,mc_cover,dep,my_lambda,my_mc=cover_for_log("c2t.log")






