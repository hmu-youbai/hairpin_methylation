# hairpin_methylation
实现cut fq1和fq2，同时可以提取hmc和mc，内测版本。

cai_club

Just for test and Processing simple data

## 安装 
ps:建议在conda虚拟环境下安装

1. 克隆代码库到本地：

   ```bash
   git clone https://github.com/hmu-youbai/hairpin_methylation.git
   ```
   
   如果网络问题无法完成克隆，可以直接下载.zip压缩包，自行解压
   
2. 进入项目目录：
 
   ```bash
   cd hairpin_methylation
   ```
3. 运行安装命令：   

   ```bash
   python setup.py build
   ```
   ```bash
   python setup.py install
   ```
   注意：确保你的系统已经安装了所需的依赖包。


## 使用 

0.cai_club内部命令

   ```bash
   run_hairpin --fq1 yourfq1.fq.gz --fq2 yourfq2.fq.gz --parallel 20 --duplication 1
   ```
   run_hairpin命令会执行全部操作：cut，提取hmc，提取mc，并生成log文件


1. 切割原始fq1和fq2，并生成3个文件：①还原到BS之前的stem 5` - 3`的fq文件，②仅切割stem的fq1，③仅切割stem的fq2

   ```bash
   hairpin_cut --fq1 yourfq1.fq.gz --fq2 yourfq2.fq.gz --parallel 20 --duplication 1
    ```
   
   如果已经自行完成去重，可以删除duplication选项
   parallel选项设置进程数，当文件过大时（超过50g），建议设置5，避免内存溢出
   默认模式为填充普通C，如果填充mC请设置--rule 1


2. ①文件比对得到的sam文件提取hmc：
 
   ```bash
   hmc_extractor --sam restored_seq.fq.sam --parallel 20
    ```
  

3. ①文件比对得到的sam文件提取mc：
 
   ```bash
   mc_extractor --sam restored_seq.fq.sam  --cutfq1 restored_seq.cut_f1.fq --parallel 20
    ```

