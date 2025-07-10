"""
GWAS数据处理模块
"""

import pandas as pd
from pathlib import Path
from .utils import CommandRunner

class DataProcessor:
    """数据处理器"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.plink_prefix_qc = None
    
    def vcf_to_plink(self):
        """VCF转换为PLINK格式"""
        vcf_path = Path(self.config.vcf_file).resolve()
        cmd = (
            f"plink --vcf {vcf_path} --make-bed --out {self.config.base_name} "
            f"--allow-extra-chr --double-id"
        )
        self.cmd_runner.run(cmd, "VCF转换为PLINK格式")
    
    def quality_control(self):
        """质量控制 | Quality control"""
        self.logger.info(
            f"使用缺失率阈值 | Using missing rate thresholds: --geno {self.config.missing_rate} (SNP), "
            f"--mind {self.config.missing_rate} (样本 | samples)"
        )
        
        # MAF质控 | MAF filtering
        cmd_maf = (
            f"plink --bfile {self.config.base_name} --maf {self.config.maf} "
            f"--make-bed --out {self.config.base_name}_maf --allow-extra-chr"
        )
        self.cmd_runner.run(cmd_maf, f"MAF质控 | MAF filtering (>={self.config.maf})")
        
        # HWE质控 | HWE filtering
        cmd_hwe = (
            f"plink --bfile {self.config.base_name}_maf --hwe {self.config.hwe_pvalue} "
            f"--make-bed --out {self.config.base_name}_hwe --allow-extra-chr"
        )
        self.cmd_runner.run(cmd_hwe, f"Hardy-Weinberg平衡检验 | Hardy-Weinberg equilibrium test (p>{self.config.hwe_pvalue})")
        
        # 缺失率质控 | Missing rate filtering
        cmd_miss = (
            f"plink --bfile {self.config.base_name}_hwe "
            f"--geno {self.config.missing_rate} --mind {self.config.missing_rate} "
            f"--make-bed --out {self.config.base_name}_qc --allow-extra-chr"
        )
        self.cmd_runner.run(cmd_miss, f"缺失率质控 | Missing rate filtering (<{self.config.missing_rate})")
        
        self.plink_prefix_qc = f"{self.config.base_name}_qc"
        self.log_qc_stats()
    
    def log_qc_stats(self):
        """记录质控统计 | Log QC statistics"""
        try:
            fam_file = self.config.output_path / f"{self.plink_prefix_qc}.fam"
            bim_file = self.config.output_path / f"{self.plink_prefix_qc}.bim"
            
            if fam_file.exists():
                sample_count = sum(1 for line in open(fam_file))
                self.logger.info(f"质控后样本数 | Samples after QC: {sample_count}")
            
            if bim_file.exists():
                snp_count = sum(1 for line in open(bim_file))
                self.logger.info(f"质控后SNP数 | SNPs after QC: {snp_count}")
                
        except Exception as e:
            self.logger.warning(f"无法读取质控统计 | Cannot read QC statistics: {e}")

class PhenotypeProcessor:
    """表型数据处理器"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def prepare_plink_pheno_file(self, plink_prefix_qc: str):
        """准备PLINK表型文件"""
        fam_path = self.config.output_path / f"{plink_prefix_qc}.fam"
        fam_df = pd.read_csv(
            fam_path, 
            sep=r'\s+', 
            header=None, 
            usecols=[0, 1], 
            names=['FID', 'IID']
        )
        
        raw_pheno_df = pd.read_csv(self.config.phenotype_file, sep='\t')
        pheno_id_col = raw_pheno_df.columns[0]
        raw_pheno_df.rename(columns={pheno_id_col: 'IID'}, inplace=True)
        
        merged_df = pd.merge(fam_df, raw_pheno_df, on='IID', how='left')
        plink_pheno_path = self.config.output_path / "plink_pheno_input.txt"
        trait_names = [col for col in raw_pheno_df.columns if col != 'IID']
        
        # 二分类表型转换 | Binary phenotype conversion
        if self.config.pheno_type == 'binary':
            self.logger.info(
                "检测到二分类表型，正在转换为PLINK标准编码 | "
                "Detected binary phenotype, converting to PLINK standard encoding "
                "(抗病/对照 0 -> 1; 发病/病例 1 -> 2 | resistant/control 0 -> 1; susceptible/case 1 -> 2)"
            )
            for trait in trait_names:
                self.logger.info(f"  - 转换性状 | Converting trait: '{trait}'")
                merged_df[trait] = merged_df[trait].map({0: 1, 1: 2})
        
        output_cols = ['FID', 'IID'] + trait_names
        merged_df[output_cols].to_csv(
            plink_pheno_path, 
            sep='\t', 
            header=False, 
            index=False, 
            na_rep='-9'
        )
        
        self.logger.info(f"已为PLINK准备好多表型文件 | Multi-phenotype file prepared for PLINK: {plink_pheno_path.name}")
        
        # 显示文件内容预览 | Show file content preview
        with open(plink_pheno_path, 'r') as f:
            preview = "".join(f.readlines()[:5]).strip()
            self.logger.info(
                f"'{plink_pheno_path.name}' 的前5行内容 | First 5 lines of '{plink_pheno_path.name}':\n---\n{preview}\n---"
            )
        
        return trait_names, plink_pheno_path
