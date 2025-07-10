"""
GWAS分析核心模块
"""

import pandas as pd
from pathlib import Path
from .utils import CommandRunner

class GWASAnalyzer:
    """GWAS分析器"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_gwas_for_trait(self, trait_name: str, trait_index: int, 
                          plink_prefix_qc: str, plink_pheno_path: Path):
        """为单个表型运行GWAS分析"""
        trait_bfile_prefix = f"{plink_prefix_qc}_{trait_name}"
        
        # 创建表型专用PLINK文件
        plink_update_cmd = (
            f"plink --bfile {plink_prefix_qc} "
            f"--pheno {plink_pheno_path.name} "
            f"--mpheno {trait_index + 1} "
            f"--make-bed --out {trait_bfile_prefix} --allow-extra-chr"
        )
        self.cmd_runner.run(
            plink_update_cmd, 
            f"为性状 '{trait_name}' 创建专用的PLINK文件集"
        )
        
        # 计算亲缘关系矩阵
        kinship_prefix = f"kinship_{trait_name}"
        kinship_cmd = f"gemma -bfile {trait_bfile_prefix} -gk 1 -o {kinship_prefix}"
        self.cmd_runner.run(kinship_cmd, "计算亲缘关系矩阵")
        
        # 处理PCA协变量
        pca_covar_arg = self._handle_pca_covariates(trait_name, trait_bfile_prefix)
        
        # 构建GWAS命令
        gwas_cmd = self._build_gwas_command(
            trait_bfile_prefix, trait_name, kinship_prefix, pca_covar_arg
        )
        
        self.cmd_runner.run(
            gwas_cmd, 
            f"GWAS分析 - {trait_name} ({self.config.model}模型)"
        )
        
        return trait_bfile_prefix
    
    def _handle_pca_covariates(self, trait_name: str, trait_bfile_prefix: str) -> str:
        """处理PCA协变量"""
        if not self.config.use_pca:
            return ""
        
        if self.config.model == "bslmm":
            self.logger.warning("BSLMM模型不支持协变量，将忽略PCA设置")
            return ""
        
        # 只有LMM模型才处理PCA
        pca_prefix = f"pca_{trait_name}"
        pca_cmd = (
            f"plink --bfile {trait_bfile_prefix} --pca {self.config.pca_components} "
            f"--out {pca_prefix} --allow-extra-chr"
        )
        self.cmd_runner.run(pca_cmd, f"PCA分析 (前{self.config.pca_components}个主成分)")
        
        pca_eigenvec_path = self.config.output_path / f"{pca_prefix}.eigenvec"
        if pca_eigenvec_path.exists():
            pca_df = pd.read_csv(pca_eigenvec_path, sep=r'\s+', header=None)
            pca_covar = pca_df.iloc[:, 2:2+self.config.pca_components]
            pca_covar_path = self.config.output_path / f"pca_covariates_{trait_name}.txt"
            pca_covar.to_csv(pca_covar_path, sep='\t', index=False, header=False)
            return f"-c {pca_covar_path.name}"
        
        return ""
    
    def _build_gwas_command(self, trait_bfile_prefix: str, trait_name: str, 
                           kinship_prefix: str, pca_covar_arg: str) -> str:
        """构建GWAS命令"""
        gwas_cmd_parts = [
            f"gemma -bfile {trait_bfile_prefix}",
            f"-k output/{kinship_prefix}.cXX.txt",
            pca_covar_arg,
            f"-o {trait_name}_gwas"
        ]
        
        if self.config.model == "lmm":
            if self.config.pheno_type == 'binary':
                self.logger.info("表型为二分类，使用LMM似然比检验 (-lmm 2)")
                gwas_cmd_parts.append("-lmm 2")
            else:
                self.logger.info("表型为数量性状，使用LMM Wald检验 (-lmm 4)")
                gwas_cmd_parts.append("-lmm 4")
        elif self.config.model == "bslmm":
            self.logger.info("使用BSLMM模型 (-bslmm 1)，它将自动适应表型类型")
            gwas_cmd_parts.append("-bslmm 1 -w 50000 -s 500000")
        
        return " ".join(filter(None, gwas_cmd_parts))
