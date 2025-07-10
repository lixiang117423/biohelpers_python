#!/usr/bin/env python3
"""
GWAS分析运行脚本 | GWAS Analysis Runner Script
这是一个简单的入口脚本，用于运行GWAS分析 | Simple entry script for running GWAS analysis

用法 | Usage:
    python run_gwas.py --vcf data.vcf --pheno pheno.txt
    
示例 | Examples:
    # 基本分析 | Basic analysis
    python run_gwas.py --vcf my_data.vcf --pheno my_pheno.txt
    
    # 使用PCA协变量的二分类表型分析 | Binary trait analysis with PCA covariates
    python run_gwas.py --vcf data.vcf --pheno pheno.txt --pheno-type binary --use-pca --pca-components 5
    
    # 使用FDR校正 | Use FDR correction
    python run_gwas.py --vcf data.vcf --pheno pheno.txt --correction fdr --sig-threshold 0.1
"""

from gwas_toolkit.main import main

if __name__ == "__main__":
    main()
