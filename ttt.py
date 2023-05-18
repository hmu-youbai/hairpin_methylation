def validate_fastq_file(fastq_file):
    """
    验证 FASTQ 文件中的序列长度和质量值长度是否一致

    参数:
        - fastq_file (str): FASTQ 文件的文件路径

    返回:
        - bool: 序列长度和质量值长度是否一致，True 表示一致，False 表示不一致
    """
    with open(fastq_file, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 4):
            # 读取每个记录的4行内容
            record_id = lines[i].strip()
            sequence = lines[i + 1].strip()
            quality = lines[i + 3].strip()

            # 检查序列长度和质量值长度是否一致
            if len(sequence) != len(quality):
                return False
    return True

# 调用示例
fastq_file = 'restored_seq.fq'  # 替换为你的 FASTQ 文件的文件路径
result = validate_fastq_file(fastq_file)
if result:
    print('序列长度和质量值长度一致')
else:
    print('序列长度和质量值长度不一致')

