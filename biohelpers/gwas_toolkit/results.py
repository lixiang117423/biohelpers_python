"""
GWAS结果处理模块
"""

import pandas as pd
from pathlib import Path

try:
    from statsmodels.stats.multitest import multipletests
except ImportError:
    multipletests = None

class ResultsProcessor:
    """结果处理器"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def filter_significant_snps(self, trait_name: str, trait_bfile_prefix: str):
        """筛选显著SNP | Filter significant SNPs"""
        result_file = self.config.output_path / f"output/{trait_name}_gwas.{'assoc.txt' if self.config.model == 'lmm' else 'param.txt'}"
        
        if not result_file.exists():
            self.logger.warning(f"未找到结果文件 | Result file not found: {result_file}")
            return
        
        try:
            gwas_results = pd.read_csv(result_file, sep='\t')
            
            if self.config.model == "lmm":
                significant_snps = self._filter_lmm_results(gwas_results, trait_bfile_prefix)
            else:  # bslmm
                significant_snps = self._filter_bslmm_results(gwas_results)
            
            if significant_snps.empty:
                self.logger.info(f"{trait_name}: 未发现显著SNP | No significant SNPs found.")
                return
            
            # 保存结果 | Save results
            sig_output = self.config.output_path / f"{trait_name}_significant_snps_{self.config.correction_method}.txt"
            significant_snps.to_csv(sig_output, sep='\t', index=False)
            self.logger.info(
                f"{trait_name}: 发现 {len(significant_snps)} 个显著SNP | Found {len(significant_snps)} significant SNPs，"
                f"已保存到 | saved to {sig_output.name}"
            )
            
        except Exception as e:
            self.logger.error(f"筛选显著SNP时发生错误 | Error occurred while filtering significant SNPs ({trait_name}): {e}")
    
    def _filter_lmm_results(self, gwas_results: pd.DataFrame, trait_bfile_prefix: str) -> pd.DataFrame:
        """筛选LMM结果"""
        p_col = 'p_lrt' if self.config.pheno_type == 'binary' else 'p_wald'
        
        if p_col not in gwas_results.columns:
            self.logger.warning(f"未找到期望的P值列 '{p_col}'")
            return pd.DataFrame()
        
        gwas_results.dropna(subset=[p_col], inplace=True)
        if gwas_results.empty:
            self.logger.info("没有有效的P值进行分析。")
            return pd.DataFrame()
        
        p_values = gwas_results[p_col]
        
        if self.config.correction_method == 'bonferroni':
            return self._bonferroni_correction(gwas_results, p_values, trait_bfile_prefix, p_col)
        elif self.config.correction_method == 'fdr':
            return self._fdr_correction(gwas_results, p_values, p_col)
        else:  # 'none'
            return self._no_correction(gwas_results, p_values, p_col)
    
    def _filter_bslmm_results(self, gwas_results: pd.DataFrame) -> pd.DataFrame:
        """筛选BSLMM结果"""
        gamma_threshold = 0.01
        self.logger.info(f"BSLMM模型: 筛选 gamma > {gamma_threshold} 的SNP (不使用P值校正)")
        
        significant_snps = gwas_results[gwas_results['gamma'] > gamma_threshold].copy()
        if not significant_snps.empty:
            significant_snps.sort_values(by='gamma', ascending=False, inplace=True)
        
        return significant_snps
    
    def _bonferroni_correction(self, gwas_results: pd.DataFrame, p_values: pd.Series, 
                              trait_bfile_prefix: str, p_col: str) -> pd.DataFrame:
        """Bonferroni校正"""
        bim_file = self.config.output_path / f"{trait_bfile_prefix}.bim"
        n_snps = sum(1 for line in open(bim_file))
        threshold = self.config.sig_threshold / n_snps if n_snps > 0 else 0
        
        self.logger.info(
            f"使用Bonferroni校正: 阈值 = {self.config.sig_threshold} / {n_snps} = {threshold:.2e}"
        )
        
        significant_snps = gwas_results[p_values < threshold].copy()
        if not significant_snps.empty:
            significant_snps.sort_values(by=p_col, inplace=True)
        
        return significant_snps
    
    def _fdr_correction(self, gwas_results: pd.DataFrame, p_values: pd.Series, p_col: str) -> pd.DataFrame:
        """FDR校正"""
        if multipletests is None:
            self.logger.error("FDR校正需要statsmodels库，请安装: pip install statsmodels")
            return pd.DataFrame()
        
        alpha_fdr = self.config.sig_threshold
        self.logger.info(f"使用FDR (Benjamini-Hochberg)校正，目标q-value = {alpha_fdr}")
        
        reject, pvals_corrected, _, _ = multipletests(
            p_values, alpha=alpha_fdr, method='fdr_bh'
        )
        gwas_results['q_value'] = pvals_corrected
        
        significant_snps = gwas_results[reject].copy()
        if not significant_snps.empty:
            significant_snps.sort_values(by=p_col, inplace=True)
        
        return significant_snps
    
    def _no_correction(self, gwas_results: pd.DataFrame, p_values: pd.Series, p_col: str) -> pd.DataFrame:
        """无校正"""
        threshold = self.config.sig_threshold
        self.logger.info(f"使用固定的显著性阈值: p < {threshold:.2e}")
        
        significant_snps = gwas_results[p_values < threshold].copy()
        if not significant_snps.empty:
            significant_snps.sort_values(by=p_col, inplace=True)
        
        return significant_snps

class SummaryGenerator:
    """总结生成器"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary(self, trait_names):
        """生成分析总结 | Generate analysis summary"""
        summary_file = self.config.output_path / "analysis_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("GEMMA GWAS分析总结 | GEMMA GWAS Analysis Summary\n" + "="*50 + "\n\n")
            
            # 输入文件信息 | Input file information
            f.write(f"输入文件 | Input Files:\n")
            f.write(f"  VCF文件 | VCF file: {Path(self.config.vcf_file).resolve()}\n")
            f.write(f"  表型文件 | Phenotype file: {Path(self.config.phenotype_file).resolve()}\n\n")
            
            # 分析参数 | Analysis parameters
            f.write(f"分析参数 | Analysis Parameters:\n")
            f.write(f"  MAF阈值 | MAF threshold: {self.config.maf}\n")
            f.write(f"  HWE p值阈值 | HWE p-value threshold: {self.config.hwe_pvalue}\n")
            f.write(f"  缺失率阈值 | Missing rate threshold: {self.config.missing_rate}\n")
            f.write(f"  模型类型 | Model type: {self.config.model}\n")
            f.write(f"  表型类型 | Phenotype type: {self.config.pheno_type}\n")
            f.write(f"  使用PCA协变量 | Use PCA covariates: {self.config.use_pca}\n")
            
            if self.config.use_pca:
                f.write(f"  PCA主成分数 | PCA components: {self.config.pca_components}\n")
            
            f.write(f"  P值校正方法 | P-value correction method: {self.config.correction_method}\n")
            f.write(f"  名义显著性水平(alpha/q-value) | Nominal significance level: {self.config.sig_threshold}\n")
            
            # 分析的表型 | Analyzed traits
            if trait_names:
                f.write(f"\n分析的表型 | Analyzed traits: {', '.join(trait_names)}\n")
            
            f.write(f"\n输出目录 | Output directory: {self.config.output_path.resolve()}\n")
        
        self.logger.info(f"分析总结已保存 | Analysis summary saved: {summary_file}")
