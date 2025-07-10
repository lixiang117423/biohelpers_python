"""
GEMMA GWAS分析工具包 | GEMMA GWAS Analysis Toolkit
功能: VCF到GWAS分析的完整流程，支持二分类和数量性状 | 
Features: Complete pipeline from VCF to GWAS analysis, supporting binary and quantitative traits
作者 | Author: Claude  
版本 | Version: v10 - 模块化重构版 | Modular refactored version
日期 | Date: 2025-07-10

使用示例 | Usage Examples:
    from gwas_toolkit import GemmaGWAS, GWASConfig
    
    # 创建分析器 | Create analyzer
    gwas = GemmaGWAS(
        vcf_file="data.vcf",
        phenotype_file="pheno.txt",
        pheno_type="binary",
        use_pca=True
    )
    
    # 运行分析 | Run analysis
    gwas.run_analysis()
"""

__version__ = "10.0.0"
__author__ = "Claude"

from .main import GemmaGWAS
from .config import GWASConfig

__all__ = ['GemmaGWAS', 'GWASConfig']
