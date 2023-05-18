import subprocess

# 输入文件名和输出文件名
input_file = "CGS-1w-1th_1.clean.fq.gz"
output_file = "CGS-1w-1th_2.clean.fq.gz"

# cmd = ["gunzip", input_file, output_file]
# subprocess.call(cmd)
temp=input_file[:-3]+".list"
# with open(temp, "w") as f:
#     subprocess.call(["ls", input_file[:-3], output_file[:-3]], stdout=f)
# cmd = ["fastuniq", "-i", temp, "-o", input_file[:-3]+".uniq", "-p", output_file[:-3]+".uniq"]
# subprocess.call(cmd)
cmd = ["rm", temp]
subprocess.call(cmd)
cmd = ["gzip", input_file[:-3],output_file[:-3], input_file[:-3]+".uniq", output_file[:-3]+".uniq"]
subprocess.call(cmd)
