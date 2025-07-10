# =============================================================================
# File: biohelpers/annovar/utils.py
# =============================================================================
"""
ANNOVAR Utility Functions
ANNOVAR工具函数

Utility functions for ANNOVAR annotation workflow.
ANNOVAR注释工作流程的工具函数。
"""

import os
import logging
import shutil
from pathlib import Path
from typing import Optional, Union, List, Dict, Any

from .exceptions import ANNOVARValidationError, ANNOVARFileNotFoundError


def setup_annovar_logging(
    log_file: Optional[str] = None,
    log_level: str = "INFO",
    console_output: bool = True
) -> logging.Logger:
    """
    Setup logging for ANNOVAR module.
    为ANNOVAR模块设置日志。
    
    Args:
        log_file (Optional[str]): Log file path / 日志文件路径
        log_level (str): Logging level / 日志级别
        console_output (bool): Enable console output / 启用控制台输出
        
    Returns:
        logging.Logger: Configured logger / 配置好的日志记录器
    """
    logger = logging.getLogger("biohelpers.annovar")
    logger.setLevel(getattr(logging, log_level.upper()))
    
    # Remove existing handlers / 移除现有处理器
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Create formatter / 创建格式化器
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Console handler / 控制台处理器
    if console_output:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    
    # File handler / 文件处理器
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


def validate_annovar_files(config: Dict[str, Any]) -> Dict[str, Path]:
    """
    Validate ANNOVAR input files and paths.
    验证ANNOVAR输入文件和路径。
    
    Args:
        config (Dict[str, Any]): Configuration dictionary / 配置字典
        
    Returns:
        Dict[str, Path]: Validated paths / 验证后的路径
        
    Raises:
        ANNOVARValidationError: If validation fails / 如果验证失败
        ANNOVARFileNotFoundError: If files don't exist / 如果文件不存在
    """
    required_files = {
        'gff3_file': 'GFF3 annotation file / GFF3注释文件',
        'genome_file': 'Genome sequence file / 基因组序列文件', 
        'vcf_file': 'VCF variant file / VCF变异文件'
    }
    
    required_dirs = {
        'annovar_path': 'ANNOVAR installation directory / ANNOVAR安装目录'
    }
    
    validated_paths = {}
    
    # Validate required files / 验证必需文件
    for key, description in required_files.items():
        if key not in config:
            raise ANNOVARValidationError(
                f"Missing required parameter: {key} ({description})\n"
                f"缺少必需参数: {key} ({description})"
            )
        
        file_path = Path(config[key]).resolve()
        if not file_path.exists():
            raise ANNOVARFileNotFoundError(str(file_path), "file")
        
        if not file_path.is_file():
            raise ANNOVARValidationError(
                f"Expected file but found directory: {file_path}\n"
                f"期望文件但找到目录: {file_path}"
            )
        
        validated_paths[key] = file_path
    
    # Validate required directories / 验证必需目录
    for key, description in required_dirs.items():
        if key not in config:
            raise ANNOVARValidationError(
                f"Missing required parameter: {key} ({description})\n"
                f"缺少必需参数: {key} ({description})"
            )
        
        dir_path = Path(config[key]).resolve()
        if not dir_path.exists():
            raise ANNOVARFileNotFoundError(str(dir_path), "directory")
        
        if not dir_path.is_dir():
            raise ANNOVARValidationError(
                f"Expected directory but found file: {dir_path}\n"
                f"期望目录但找到文件: {dir_path}"
            )
        
        validated_paths[key] = dir_path
    
    return validated_paths


def check_annovar_dependencies() -> Dict[str, bool]:
    """
    Check if required ANNOVAR dependencies are available.
    检查必需的ANNOVAR依赖是否可用。
    
    Returns:
        Dict[str, bool]: Availability status of each tool / 各工具的可用状态
    """
    dependencies = {
        'gff3ToGenePred': 'UCSC genome browser tool for GFF3 conversion',
        'bcftools': 'BCFtools for VCF processing (optional)',
        'perl': 'Perl interpreter for ANNOVAR scripts'
    }
    
    status = {}
    
    for tool, description in dependencies.items():
        status[tool] = shutil.which(tool) is not None
    
    return status


def normalize_build_version(build_ver: str) -> str:
    """
    Normalize and validate build version string.
    标准化并验证构建版本字符串。
    
    Args:
        build_ver (str): Raw build version / 原始构建版本
        
    Returns:
        str: Normalized build version / 标准化的构建版本
        
    Raises:
        ANNOVARValidationError: If build version is invalid / 如果构建版本无效
    """
    if not build_ver:
        raise ANNOVARValidationError(
            "Build version cannot be empty\n"
            "构建版本不能为空"
        )
    
    # Remove path separators / 移除路径分隔符
    clean_build_ver = os.path.basename(build_ver.strip())
    
    # Check for invalid characters / 检查无效字符
    invalid_chars = ['/', '\\', ' ', '\t', '\n']
    if any(char in clean_build_ver for char in invalid_chars):
        raise ANNOVARValidationError(
            f"Build version contains invalid characters: {build_ver}\n"
            f"构建版本包含无效字符: {build_ver}"
        )
    
    return clean_build_ver


def ensure_output_directory(output_dir: Union[str, Path]) -> Path:
    """
    Ensure output directory exists.
    确保输出目录存在。
    
    Args:
        output_dir (Union[str, Path]): Output directory path / 输出目录路径
        
    Returns:
        Path: Resolved output directory path / 解析后的输出目录路径
    """
    dir_path = Path(output_dir).resolve()
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def get_file_basename(file_path: Union[str, Path], extensions: Optional[List[str]] = None) -> str:
    """
    Get file basename without specified extensions.
    获取不包含指定扩展名的文件基础名称。
    
    Args:
        file_path (Union[str, Path]): File path / 文件路径
        extensions (Optional[List[str]]): Extensions to remove / 要移除的扩展名
        
    Returns:
        str: File basename / 文件基础名称
    """
    if extensions is None:
        extensions = ['.vcf', '.gz', '.txt', '.gff3', '.fa', '.fasta']
    
    basename = Path(file_path).name
    
    for ext in extensions:
        if basename.endswith(ext):
            basename = basename[:-len(ext)]
    
    return basename
