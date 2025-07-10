# =============================================================================
# File: biohelpers/annovar/cli.py
# =============================================================================
"""
Command Line Interface for ANNOVAR Module
ANNOVAR模块的命令行接口

Command line interface for the biohelpers ANNOVAR annotation module.
biohelpers ANNOVAR注释模块的命令行接口。
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, Any, Optional

from .annotator import ANNOVARAnnotator
from .exceptions import ANNOVARError
from .utils import setup_annovar_logging, check_annovar_dependencies


def create_annovar_parser() -> argparse.ArgumentParser:
    """
    Create argument parser for ANNOVAR CLI.
    为ANNOVAR CLI创建参数解析器。
    
    Returns:
        argparse.ArgumentParser: Configured argument parser / 配置好的参数解析器
    """
    parser = argparse.ArgumentParser(
        prog='biohelpers-annovar',
        description=(
            'ANNOVAR VCF Annotation Module for biohelpers\n'
            'biohelpers的ANNOVAR VCF注释模块\n\n'
            'Automated ANNOVAR annotation pipeline for VCF files.\n'
            '用于VCF文件的自动化ANNOVAR注释流水线。'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples / 示例:\n'
            '  # Full pipeline / 完整流水线:\n'
            '  biohelpers-annovar --gff3 genome.gff3 --genome genome.fa --vcf variants.vcf --build-ver OV\n\n'
            '  # Single step / 单步执行:\n'
            '  biohelpers-annovar --gff3 genome.gff3 --genome genome.fa --vcf variants.vcf --build-ver OV --step 1\n\n'
            '  # With custom output directory / 自定义输出目录:\n'
            '  biohelpers-annovar --gff3 genome.gff3 --genome genome.fa --vcf variants.vcf --build-ver OV --output-dir /path/to/output\n'
        )
    )
    
    # Required arguments group / 必需参数组
    required = parser.add_argument_group('Required Arguments / 必需参数')
    
    required.add_argument(
        '--gff3',
        required=True,
        metavar='FILE',
        help=(
            'GFF3 annotation file path\n'
            'GFF3注释文件路径'
        )
    )
    
    required.add_argument(
        '--genome',
        required=True,
        metavar='FILE',
        help=(
            'Genome sequence file path (FASTA format)\n'
            '基因组序列文件路径（FASTA格式）'
        )
    )
    
    required.add_argument(
        '--vcf',
        required=True,
        metavar='FILE',
        help=(
            'VCF variant file path\n'
            'VCF变异文件路径'
        )
    )
    
    required.add_argument(
        '--build-ver',
        required=True,
        metavar='VERSION',
        help=(
            'Genome build version identifier (e.g., hg38, OV, KY131)\n'
            '基因组构建版本标识符（例如：hg38, OV, KY131）\n'
            'Must not contain path separators\n'
            '不得包含路径分隔符'
        )
    )
    
    # Optional arguments group / 可选参数组
    optional = parser.add_argument_group('Optional Arguments / 可选参数')
    
    optional.add_argument(
        '--annovar-path',
        metavar='DIR',
        default='/share/org/YZWL/yzwl_lixg/software/annovar/annovar',
        help=(
            'ANNOVAR installation directory path\n'
            'ANNOVAR安装目录路径\n'
            '(default: /share/org/YZWL/yzwl_lixg/software/annovar/annovar)\n'
            '（默认值：/share/org/YZWL/yzwl_lixg/software/annovar/annovar）'
        )
    )
    
    optional.add_argument(
        '--database-path',
        metavar='DIR',
        default='./database',
        help=(
            'ANNOVAR database directory path\n'
            'ANNOVAR数据库目录路径\n'
            '(default: ./database)\n'
            '（默认值：./database）'
        )
    )
    
    optional.add_argument(
        '--output-dir',
        metavar='DIR',
        default='./annovar_output',
        help=(
            'Output directory path\n'
            '输出目录路径\n'
            '(default: ./annovar_output)\n'
            '（默认值：./annovar_output）'
        )
    )
    
    # Step control group / 步骤控制组
    step_control = parser.add_argument_group('Step Control / 步骤控制')
    
    step_control.add_argument(
        '--step',
        type=int,
        choices=[1, 2, 3, 4],
        metavar='N',
        help=(
            'Execute only specified step (1-4)\n'
            '只执行指定的步骤（1-4）\n'
            '1: GFF3 to GenePred conversion / GFF3到GenePred转换\n'
            '2: Extract transcript sequences / 提取转录本序列\n'
            '3: VCF processing and conversion / VCF处理和转换\n'
            '4: Variant annotation / 变异注释'
        )
    )
    
    # VCF processing group / VCF处理组
    vcf_group = parser.add_argument_group('VCF Processing Options / VCF处理选项')
    
    vcf_group.add_argument(
        '--qual-threshold',
        type=int,
        default=20,
        metavar='N',
        help=(
            'VCF quality filtering threshold\n'
            'VCF质量过滤阈值\n'
            '(default: 20, only used when filtering is enabled)\n'
            '（默认值：20，仅在启用过滤时使用）'
        )
    )
    
    vcf_filter_group = vcf_group.add_mutually_exclusive_group()
    
    vcf_filter_group.add_argument(
        '--skip-vcf-filter',
        action='store_true',
        default=True,
        help=(
            'Skip VCF quality filtering (default)\n'
            '跳过VCF质量过滤（默认）'
        )
    )
    
    vcf_filter_group.add_argument(
        '--enable-vcf-filter',
        action='store_true',
        help=(
            'Enable VCF quality filtering using bcftools\n'
            '启用使用bcftools的VCF质量过滤'
        )
    )
    
    # GFF3 processing group / GFF3处理组
    gff_group = parser.add_argument_group('GFF3 Processing Options / GFF3处理选项')
    
    gff_group.add_argument(
        '--skip-gff-fix',
        action='store_true',
        help=(
            'Skip automatic GFF3 file fixes (header, CDS phase, etc.)\n'
            '跳过自动GFF3文件修复（头部、CDS phase等）'
        )
    )
    
    # Logging group / 日志组
    logging_group = parser.add_argument_group('Logging Options / 日志选项')
    
    logging_group.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help=(
            'Set logging level\n'
            '设置日志级别\n'
            '(default: INFO)\n'
            '（默认值：INFO）'
        )
    )
    
    logging_group.add_argument(
        '--log-file',
        metavar='FILE',
        help=(
            'Write logs to specified file (in addition to console)\n'
            '将日志写入指定文件（除控制台外）'
        )
    )
    
    # Utility options / 工具选项
    utility_group = parser.add_argument_group('Utility Options / 工具选项')
    
    utility_group.add_argument(
        '--check-deps',
        action='store_true',
        help=(
            'Check for required external dependencies and exit\n'
            '检查必需的外部依赖项并退出'
        )
    )
    
    utility_group.add_argument(
        '--version',
        action='version',
        version='biohelpers.annovar 1.0.0'
    )
    
    return parser


def check_dependencies_and_exit() -> None:
    """
    Check external dependencies and exit.
    检查外部依赖项并退出。
    """
    print("Checking ANNOVAR dependencies...")
    print("检查ANNOVAR依赖项...")
    print()
    
    deps = check_annovar_dependencies()
    
    all_available = True
    for tool, available in deps.items():
        status = "✓ Available" if available else "✗ Missing"
        status_cn = "✓ 可用" if available else "✗ 缺失"
        print(f"{tool}: {status} / {status_cn}")
        
        if not available:
            all_available = False
    
    print()
    if all_available:
        print("All dependencies are available! / 所有依赖项都可用!")
        sys.exit(0)
    else:
        print("Some dependencies are missing. Please install them before running ANNOVAR.")
        print("缺少一些依赖项。请在运行ANNOVAR之前安装它们。")
        sys.exit(1)


def build_config_from_args(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Build configuration dictionary from parsed arguments.
    从解析的参数构建配置字典。
    
    Args:
        args (argparse.Namespace): Parsed command line arguments / 解析的命令行参数
        
    Returns:
        Dict[str, Any]: Configuration dictionary / 配置字典
    """
    config = {
        'gff3_file': args.gff3,
        'genome_file': args.genome,
        'vcf_file': args.vcf,
        'build_ver': args.build_ver,
        'annovar_path': args.annovar_path,
        'database_path': args.database_path,
        'output_dir': args.output_dir,
        'qual_threshold': args.qual_threshold,
        'skip_vcf_filter': args.skip_vcf_filter and not args.enable_vcf_filter,
        'skip_gff_fix': args.skip_gff_fix
    }
    
    return config


def main() -> int:
    """
    Main CLI entry point for ANNOVAR module.
    ANNOVAR模块的主要CLI入口点。
    
    Returns:
        int: Exit code / 退出代码
    """
    parser = create_annovar_parser()
    args = parser.parse_args()
    
    # Handle dependency check / 处理依赖检查
    if args.check_deps:
        check_dependencies_and_exit()
    
    # Setup logging / 设置日志
    logger = setup_annovar_logging(
        log_file=args.log_file,
        log_level=args.log_level,
        console_output=True
    )
    
    # Build configuration / 构建配置
    config = build_config_from_args(args)
    
    try:
        # Initialize annotator / 初始化注释器
        logger.info("Initializing ANNOVAR annotator...")
        logger.info("正在初始化ANNOVAR注释器...")
        
        annotator = ANNOVARAnnotator(config, logger)
        
        # Execute annotation / 执行注释
        if args.step:
            # Run single step / 运行单个步骤
            logger.info(f"Running single step: {args.step}")
            logger.info(f"运行单个步骤: {args.step}")
            
            success = annotator.run_single_step(args.step)
            
            if success:
                print(f"\nStep {args.step} completed successfully!")
                print(f"步骤 {args.step} 成功完成!")
                return 0
            else:
                print(f"\nStep {args.step} failed!")
                print(f"步骤 {args.step} 失败!")
                return 1
        else:
            # Run full pipeline / 运行完整流水线
            logger.info("Running full annotation pipeline...")
            logger.info("运行完整注释流水线...")
            
            success = annotator.run_full_pipeline()
            
            if success:
                output_files = annotator.get_output_files()
                print(f"\nAnnotation pipeline completed successfully!")
                print(f"注释流水线成功完成!")
                print(f"Generated {len(output_files)} output files.")
                print(f"生成了 {len(output_files)} 个输出文件。")
                return 0
            else:
                print(f"\nAnnotation pipeline failed!")
                print(f"注释流水线失败!")
                return 1
    
    except ANNOVARError as e:
        logger.error(f"ANNOVAR error: {e}")
        print(f"\nError: {e}", file=sys.stderr)
        return 1
    
    except KeyboardInterrupt:
        logger.info("Process interrupted by user")
        logger.info("进程被用户中断")
        print("\nProcess interrupted by user / 进程被用户中断")
        return 130
    
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        print(f"\nUnexpected error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
