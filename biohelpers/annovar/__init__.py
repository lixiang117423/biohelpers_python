# =============================================================================
# File: biohelpers/annovar/__init__.py
# =============================================================================
"""
ANNOVAR Annotation Module for biohelpers
biohelpers的ANNOVAR注释模块

ANNOVAR VCF annotation workflow automation module.
ANNOVAR VCF注释工作流程自动化模块。

This module provides automated ANNOVAR annotation pipeline including:
该模块提供自动化的ANNOVAR注释流水线，包括：
- GFF3 to GenePred conversion / GFF3到GenePred转换
- Transcript sequence extraction / 转录本序列提取  
- VCF processing and conversion / VCF处理和转换
- Variant annotation / 变异注释

Example usage / 使用示例:
    from biohelpers.annovar import ANNOVARAnnotator
    
    config = {
        'gff3_file': 'genome.gff3',
        'genome_file': 'genome.fasta',
        'vcf_file': 'variants.vcf',
        'build_ver': 'OV'
    }
    
    annotator = ANNOVARAnnotator(config)
    success = annotator.run_full_pipeline()
"""

from .annotator import ANNOVARAnnotator
from .exceptions import (
    ANNOVARError,
    ANNOVARValidationError, 
    ANNOVARFileNotFoundError,
    ANNOVARCommandError
)
from .utils import (
    setup_annovar_logging,
    validate_annovar_files,
    check_annovar_dependencies
)

__version__ = "1.0.0"
__author__ = "biohelpers team"

__all__ = [
    "ANNOVARAnnotator",
    "ANNOVARError",
    "ANNOVARValidationError",
    "ANNOVARFileNotFoundError", 
    "ANNOVARCommandError",
    "setup_annovar_logging",
    "validate_annovar_files",
    "check_annovar_dependencies"
]
