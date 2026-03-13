#!/usr/bin/env python3
"""
FastQC 质量检查运行脚本
用法：python fastqc_runner.py <fastq文件>
"""

import os
import sys
import subprocess

def run_fastqc(fastq_file):
    """运行 FastQC 对 FASTQ 文件进行质量检查"""
    
    # 检查文件是否存在
    if not os.path.exists(fastq_file):
        print(f"错误：文件 {fastq_file} 不存在")
        return False
    
    # 检查 FastQC 是否已安装
    try:
        subprocess.run(['fastqc', '--version'], check=True, capture_output=True)
    except:
        print("错误：FastQC 未安装，请先安装：sudo apt install fastqc")
        return False
    
    # 运行 FastQC
    print(f"正在对 {fastq_file} 运行 FastQC...")
    try:
        subprocess.run(['fastqc', fastq_file], check=True)
        print("FastQC 运行完成！")
        return True
    except subprocess.CalledProcessError as e:
        print(f"运行失败：{e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法：python fastqc_runner.py <fastq文件>")
        sys.exit(1)
    
    fastq_file = sys.argv[1]
    success = run_fastqc(fastq_file)
    
    if not success:
        sys.exit(1)