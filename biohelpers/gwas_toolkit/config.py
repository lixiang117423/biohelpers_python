"""
GWAS分析配置管理模块 | GWAS Analysis Configuration Management Module
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class GWASConfig:
    """GWAS分析配置类 | GWAS Analysis Configuration Class"""
    
    # 输入文件
    vcf_file: str
    phenotype_file: str
    output_dir: str = "gwas_output"
    
    # 质量控制参数
    maf: float = 0.05
    hwe_pvalue: float = 1e-6
    missing_rate: float = 0.1
    
    # 分析参数
    sig_threshold: float = 1e-5
    use_pca: bool = False
    pca_components: int = 3
    model: str = "lmm"  # lmm 或 bslmm
    correction_method: str = "none"  # bonferroni, fdr, none
    pheno_type: str = "quantitative"  # binary 或 quantitative
    
    # 内部属性
    base_name: str = "gwas_data"
    
    def __post_init__(self):
        """初始化后处理"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        if not Path(self.vcf_file).exists():
            errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")
        
        if not Path(self.phenotype_file).exists():
            errors.append(f"表型文件不存在 | Phenotype file does not exist: {self.phenotype_file}")
        
        if self.model not in ["lmm", "bslmm"]:
            errors.append(f"不支持的模型类型 | Unsupported model type: {self.model}")
        
        if self.pheno_type not in ["binary", "quantitative"]:
            errors.append(f"不支持的表型类型 | Unsupported phenotype type: {self.pheno_type}")
        
        if self.correction_method not in ["bonferroni", "fdr", "none"]:
            errors.append(f"不支持的校正方法 | Unsupported correction method: {self.correction_method}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
