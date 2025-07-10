"""
GWAS分析主程序模块 | GWAS Analysis Main Module
"""

import argparse
import sys
from .config import GWASConfig
from .utils import GWASLogger, CommandRunner, check_dependencies
from .data_processing import DataProcessor, PhenotypeProcessor
from .gwas_analysis import GWASAnalyzer
from .results import ResultsProcessor, SummaryGenerator

class GemmaGWAS:
    """GWAS分析主类 | Main GWAS Analysis Class"""
    
    def __init__(self, **kwargs):
        # 检查依赖 | Check dependencies
        check_dependencies()
        
        # 初始化配置 | Initialize configuration
        self.config = GWASConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = GWASLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.data_processor = DataProcessor(self.config, self.logger, self.cmd_runner)
        self.pheno_processor = PhenotypeProcessor(self.config, self.logger)
        self.gwas_analyzer = GWASAnalyzer(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的GWAS分析流程 | Run complete GWAS analysis pipeline"""
        try:
            self.logger.info("="*20 + " 开始GEMMA GWAS分析流程 | Starting GEMMA GWAS Analysis Pipeline " + "="*20)
            
            # 步骤1: 数据预处理 | Step 1: Data preprocessing
            self.logger.info("步骤 1/3: VCF转换与质量控制 | Step 1/3: VCF conversion and quality control")
            self.data_processor.vcf_to_plink()
            self.data_processor.quality_control()
            
            # 步骤2: 表型准备 | Step 2: Phenotype preparation
            self.logger.info("\n步骤 2/3: 准备多表型文件 | Step 2/3: Preparing multi-phenotype files")
            trait_names, plink_pheno_path = self.pheno_processor.prepare_plink_pheno_file(
                self.data_processor.plink_prefix_qc
            )
            
            # 步骤3: GWAS分析 | Step 3: GWAS analysis
            self.logger.info(f"\n步骤 3/3: 循环为 {len(trait_names)} 个表型运行GWAS | Step 3/3: Running GWAS for {len(trait_names)} traits")
            for i, trait_name in enumerate(trait_names):
                self.logger.info(f"\n--- 开始处理表型 {i+1}/{len(trait_names)}: {trait_name} | Processing trait {i+1}/{len(trait_names)}: {trait_name} ---")
                
                # 运行GWAS | Run GWAS
                trait_bfile_prefix = self.gwas_analyzer.run_gwas_for_trait(
                    trait_name, i, self.data_processor.plink_prefix_qc, plink_pheno_path
                )
                
                # 筛选显著结果 | Filter significant results
                self.results_processor.filter_significant_snps(trait_name, trait_bfile_prefix)
            
            # 生成总结 | Generate summary
            self.logger.info("\n" + "="*20 + " GWAS分析全部完成 | GWAS Analysis Completed " + "="*20)
            self.summary_generator.generate_summary(trait_names)
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="GEMMA GWAS分析简化脚本 (模块化版本) | GEMMA GWAS Analysis Simplified Script (Modular Version)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 必需参数 | Required arguments
    parser.add_argument("--vcf", required=True, 
                       help="输入VCF文件 | Input VCF file")
    parser.add_argument("--pheno", required=True, 
                       help="表型文件 (tab分隔,第一列为'sample',有表头) | Phenotype file (tab-separated, first column 'sample', with header)")
    
    # 可选参数 | Optional arguments
    parser.add_argument("--output", default="gwas_output", 
                       help="输出目录 | Output directory")
    parser.add_argument("--maf", type=float, default=0.05, 
                       help="MAF阈值 | Minor allele frequency threshold")
    parser.add_argument("--hwe", type=float, default=0.001, 
                       help="HWE p值阈值 | Hardy-Weinberg equilibrium p-value threshold")
    parser.add_argument("--missing", type=float, default=0.9, 
                       help="SNP和样本的最大允许缺失率。如果数据质量较差，可能需要放宽此值(例如 0.9) | "
                            "Maximum allowed missing rate for SNPs and samples. For poor quality data, may need to relax (e.g., 0.9)")
    parser.add_argument("--sig-threshold", type=float, default=0.05,
                       help="用于P值校正的名义alpha水平(Bonferroni)或目标q值(FDR)，或固定的P值阈值 | "
                            "Nominal alpha level for p-value correction (Bonferroni) or target q-value (FDR), or fixed p-value threshold")
    parser.add_argument("--use-pca", action="store_true", 
                       help="使用PCA结果作为协变量以校正群体结构 | Use PCA results as covariates to correct for population structure")
    parser.add_argument("--pca-components", type=int, default=3, 
                       help="PCA主成分数量 | Number of PCA components")
    parser.add_argument("--model", choices=["lmm", "bslmm"], default="lmm", 
                       help="GWAS模型类型 | GWAS model type")
    parser.add_argument("--correction", choices=["bonferroni", "fdr", "none"], default="none",
                       help="P值校正方法 | P-value correction method: "
                            "'bonferroni': 严格的FWER控制 | strict FWER control; "
                            "'fdr': 控制错误发现率 | false discovery rate control; "
                            "'none': 使用固定的显著性阈值 | use fixed significance threshold")
    parser.add_argument("--pheno-type", choices=["binary", "quantitative"], default="quantitative",
                       help="表型的数据类型 | Phenotype data type: "
                            "'binary': 用于病例/对照(0/1->1/2)研究 | for case/control (0/1->1/2) studies; "
                            "'quantitative': 用于连续值性状 | for continuous traits")
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    gwas_analyzer = GemmaGWAS(**vars(args))
    gwas_analyzer.run_analysis()

if __name__ == "__main__":
    main()
