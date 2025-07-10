# =============================================================================
# File: biohelpers/annovar/annotator.py
# =============================================================================
"""
ANNOVAR Annotator Core Module
ANNOVAR注释器核心模块

Main ANNOVAR annotation workflow implementation.
主要的ANNOVAR注释工作流程实现。
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
import logging

from .exceptions import (
    ANNOVARError,
    ANNOVARValidationError,
    ANNOVARFileNotFoundError,
    ANNOVARCommandError
)
from .utils import (
    setup_annovar_logging,
    validate_annovar_files,
    normalize_build_version,
    ensure_output_directory,
    get_file_basename
)


class ANNOVARAnnotator:
    """
    ANNOVAR Annotation Workflow Manager
    ANNOVAR注释工作流程管理器
    
    Manages the complete ANNOVAR annotation pipeline including:
    管理完整的ANNOVAR注释流水线，包括：
    1. GFF3 to GenePred conversion / GFF3到GenePred转换
    2. Transcript sequence extraction / 转录本序列提取
    3. VCF processing and format conversion / VCF处理和格式转换
    4. Variant annotation / 变异注释
    
    Example:
        config = {
            'gff3_file': 'genome.gff3',
            'genome_file': 'genome.fasta', 
            'vcf_file': 'variants.vcf',
            'build_ver': 'OV',
            'annovar_path': '/path/to/annovar'
        }
        
        annotator = ANNOVARAnnotator(config)
        success = annotator.run_full_pipeline()
    """
    
    def __init__(self, config: Dict[str, Any], logger: Optional[logging.Logger] = None):
        """
        Initialize ANNOVAR annotator.
        初始化ANNOVAR注释器。
        
        Args:
            config (Dict[str, Any]): Configuration parameters / 配置参数
            logger (Optional[logging.Logger]): Custom logger / 自定义日志记录器
            
        Raises:
            ANNOVARValidationError: If configuration is invalid / 如果配置无效
            ANNOVARFileNotFoundError: If required files don't exist / 如果必需文件不存在
        """
        self.config = config.copy()
        self.logger = logger or setup_annovar_logging()
        self._intermediate_files = {}
        self._output_files = []
        
        # Validate and normalize configuration / 验证并标准化配置
        self._validate_and_setup_config()
        
        self.logger.info(
            "ANNOVAR Annotator initialized successfully.\n"
            "ANNOVAR注释器初始化成功。"
        )
    
    def _validate_and_setup_config(self) -> None:
        """
        Validate and setup configuration parameters.
        验证并设置配置参数。
        """
        # Validate required files and paths / 验证必需文件和路径
        validated_paths = validate_annovar_files(self.config)
        self.config.update(validated_paths)
        
        # Validate and normalize build version / 验证并标准化构建版本
        if 'build_ver' not in self.config:
            raise ANNOVARValidationError(
                "Missing required parameter: build_ver (Genome build version)\n"
                "缺少必需参数: build_ver (基因组构建版本)"
            )
        
        self.config['build_ver'] = normalize_build_version(self.config['build_ver'])
        
        # Setup output directory / 设置输出目录
        output_dir = self.config.get('output_dir', './annovar_output')
        self.config['output_dir'] = ensure_output_directory(output_dir)
        
        # Setup database path with default / 设置数据库路径（含默认值）
        if 'database_path' not in self.config:
            self.config['database_path'] = Path('./database').resolve()
        else:
            self.config['database_path'] = Path(self.config['database_path']).resolve()
        
        # Set default optional parameters / 设置默认可选参数
        self.config.setdefault('qual_threshold', 20)
        self.config.setdefault('skip_vcf_filter', True)
        self.config.setdefault('skip_gff_fix', False)
        
        # Log configuration summary / 记录配置摘要
        self.logger.info(f"Output directory / 输出目录: {self.config['output_dir']}")
        self.logger.info(f"Build version / 构建版本: {self.config['build_ver']}")
        self.logger.info(f"ANNOVAR path / ANNOVAR路径: {self.config['annovar_path']}")
    
    def _execute_command(self, command: str, description: str, 
                        cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
        """
        Execute system command with comprehensive error handling.
        执行系统命令并进行全面的错误处理。
        
        Args:
            command (str): Command to execute / 要执行的命令
            description (str): Command description for logging / 用于日志的命令描述
            cwd (Optional[Path]): Working directory / 工作目录
            
        Returns:
            subprocess.CompletedProcess: Command execution result / 命令执行结果
            
        Raises:
            ANNOVARCommandError: If command execution fails / 如果命令执行失败
        """
        self.logger.info(f"Executing: {description}")
        self.logger.info(f"正在执行: {description}")
        self.logger.debug(f"Command: {command}")
        
        try:
            result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=cwd
            )
            
            self.logger.info(f"Command completed successfully: {description}")
            self.logger.info(f"命令执行成功: {description}")
            
            if result.stdout:
                self.logger.debug(f"STDOUT: {result.stdout}")
            
            return result
            
        except subprocess.CalledProcessError as e:
            error_msg = (
                f"Command execution failed: {description}\n"
                f"命令执行失败: {description}\n"
                f"Command: {command}\n"
                f"命令: {command}\n"
                f"Return code: {e.returncode}\n"
                f"返回代码: {e.returncode}"
            )
            
            if e.stderr:
                error_msg += f"\nStderr: {e.stderr}\n标准错误: {e.stderr}"
            
            self.logger.error(error_msg)
            raise ANNOVARCommandError(
                error_msg,
                command=command,
                return_code=e.returncode,
                stderr=e.stderr
            )
    
    def step1_convert_gff3_to_genepred(self) -> bool:
        """
        Step 1: Convert GFF3 annotation to GenePred format.
        步骤1: 将GFF3注释转换为GenePred格式。
        
        Returns:
            bool: Success status / 成功状态
        """
        gff3_file = self.config['gff3_file']
        build_ver = self.config['build_ver']
        output_dir = self.config['output_dir']
        output_file = output_dir / f"{build_ver}_refGene.txt"
        
        self.logger.info(f"Converting GFF3 to GenePred format...")
        self.logger.info(f"正在转换GFF3到GenePred格式...")
        self.logger.info(f"Input GFF3: {gff3_file}")
        self.logger.info(f"输入GFF3: {gff3_file}")
        self.logger.info(f"Output GenePred: {output_file}")
        self.logger.info(f"输出GenePred: {output_file}")
        
        # Fix GFF3 file if needed / 如需要则修复GFF3文件
        if not self.config.get('skip_gff_fix', False):
            self._fix_gff3_issues(gff3_file)
        
        # Convert GFF3 to GenePred / 转换GFF3到GenePred
        command = f"gff3ToGenePred {gff3_file} {output_file}"
        
        try:
            self._execute_command(
                command, 
                "GFF3 to GenePred conversion / GFF3到GenePred转换"
            )
            self._intermediate_files['genepred_file'] = output_file
            return True
        except ANNOVARCommandError:
            return False
    
    def _fix_gff3_issues(self, gff3_file: Path) -> None:
        """
        Fix common issues in GFF3 files.
        修复GFF3文件中的常见问题。
        
        Args:
            gff3_file (Path): Path to GFF3 file / GFF3文件路径
        """
        self.logger.info("Checking and fixing GFF3 file issues...")
        self.logger.info("检查并修复GFF3文件问题...")
        
        # Check and add GFF3 header if missing / 检查并添加缺失的GFF3头部
        self._ensure_gff3_header(gff3_file)
        
        # Fix CDS phase issues / 修复CDS phase问题
        self._fix_cds_phase_issues(gff3_file)
    
    def _ensure_gff3_header(self, gff3_file: Path) -> None:
        """Ensure GFF3 file has proper header / 确保GFF3文件有正确的头部"""
        with open(gff3_file, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()
        
        if not first_line.startswith('##gff-version'):
            self.logger.warning("Adding missing GFF3 version header...")
            self.logger.warning("添加缺失的GFF3版本头部...")
            
            temp_file = gff3_file.with_suffix('.tmp')
            with open(gff3_file, 'r', encoding='utf-8') as original, \
                 open(temp_file, 'w', encoding='utf-8') as temp:
                temp.write('##gff-version 3\n')
                temp.write(original.read())
            
            temp_file.replace(gff3_file)
            self.logger.info("GFF3 version header added successfully.")
            self.logger.info("GFF3版本头部添加成功。")
    
    def _fix_cds_phase_issues(self, gff3_file: Path) -> None:
        """Fix CDS phase issues in GFF3 file / 修复GFF3文件中的CDS phase问题"""
        self.logger.info("Fixing CDS phase issues in GFF3 file...")
        self.logger.info("修复GFF3文件中的CDS phase问题...")
        
        temp_file = gff3_file.with_suffix('.phase_fixed')
        lines_fixed = 0
        
        with open(gff3_file, 'r', encoding='utf-8') as input_file, \
             open(temp_file, 'w', encoding='utf-8') as output_file:
            
            for line_num, line in enumerate(input_file, 1):
                line = line.rstrip('\n\r')
                
                # Skip comments and empty lines / 跳过注释和空行
                if line.startswith('#') or not line.strip():
                    output_file.write(line + '\n')
                    continue
                
                fields = line.split('\t')
                
                # Process CDS features / 处理CDS特征
                if len(fields) >= 3 and fields[2].upper() == 'CDS':
                    fixed_fields = self._fix_cds_line(fields, line_num)
                    if fixed_fields != fields:
                        lines_fixed += 1
                    output_file.write('\t'.join(fixed_fields) + '\n')
                else:
                    output_file.write(line + '\n')
        
        if lines_fixed > 0:
            # Backup and replace original file / 备份并替换原文件
            backup_file = gff3_file.with_suffix('.backup')
            gff3_file.rename(backup_file)
            temp_file.rename(gff3_file)
            
            self.logger.info(f"Fixed {lines_fixed} CDS phase issues.")
            self.logger.info(f"修复了 {lines_fixed} 个CDS phase问题。")
            self.logger.info(f"Original file backed up to: {backup_file}")
            self.logger.info(f"原文件已备份到: {backup_file}")
        else:
            temp_file.unlink()
            self.logger.info("No CDS phase issues found.")
            self.logger.info("未发现CDS phase问题。")
    
    def _fix_cds_line(self, fields: List[str], line_num: int) -> List[str]:
        """Fix CDS line fields / 修复CDS行字段"""
        # Ensure we have at least 9 fields / 确保至少有9个字段
        while len(fields) < 9:
            fields.append('.')
        
        # Fix phase field (index 7) / 修复phase字段（索引7）
        if len(fields) > 7:
            phase = fields[7]
            if phase not in ['0', '1', '2']:
                fields[7] = '0'
                self.logger.debug(f"Line {line_num}: Fixed invalid phase '{phase}' to '0'")
                self.logger.debug(f"行 {line_num}: 修复无效phase '{phase}' 为 '0'")
        
        return fields
    
    def step2_extract_transcript_sequences(self) -> bool:
        """
        Step 2: Extract transcript sequences using ANNOVAR.
        步骤2: 使用ANNOVAR提取转录本序列。
        
        Returns:
            bool: Success status / 成功状态
        """
        if 'genepred_file' not in self._intermediate_files:
            raise ANNOVARError(
                "GenePred file not available. Please run step 1 first.\n"
                "GenePred文件不可用。请先运行步骤1。"
            )
        
        genome_file = self.config['genome_file']
        genepred_file = self._intermediate_files['genepred_file']
        build_ver = self.config['build_ver']
        output_dir = self.config['output_dir']
        output_file = output_dir / f"{build_ver}_refGeneMrna.fa"
        annovar_path = self.config['annovar_path']
        
        self.logger.info("Extracting transcript sequences...")
        self.logger.info("正在提取转录本序列...")
        
        command = (
            f"perl {annovar_path}/retrieve_seq_from_fasta.pl "
            f"--format refGene --seqfile {genome_file} "
            f"{genepred_file} -outfile {output_file}"
        )
        
        try:
            self._execute_command(
                command,
                "Transcript sequence extraction / 转录本序列提取"
            )
            self._intermediate_files['mrna_file'] = output_file
            return True
        except ANNOVARCommandError:
            return False
    
    def step3_process_and_convert_vcf(self) -> bool:
        """
        Step 3: Process and convert VCF file to ANNOVAR format.
        步骤3: 处理并转换VCF文件为ANNOVAR格式。
        
        Returns:
            bool: Success status / 成功状态
        """
        vcf_file = self.config['vcf_file']
        skip_filter = self.config.get('skip_vcf_filter', True)
        output_dir = self.config['output_dir']
        annovar_path = self.config['annovar_path']
        
        vcf_basename = get_file_basename(vcf_file)
        
        # Step 3a: Filter VCF if requested / 步骤3a: 如果需要则过滤VCF
        if skip_filter:
            self.logger.info("Skipping VCF filtering (using original file)")
            self.logger.info("跳过VCF过滤（使用原始文件）")
            filtered_vcf = vcf_file
        else:
            filtered_vcf = self._filter_vcf_file(vcf_file, vcf_basename)
            if not filtered_vcf:
                return False
        
        # Step 3b: Convert to ANNOVAR format / 步骤3b: 转换为ANNOVAR格式
        return self._convert_vcf_to_annovar_format(filtered_vcf, vcf_basename)
    
    def _filter_vcf_file(self, vcf_file: Path, basename: str) -> Optional[Path]:
        """Filter VCF file by quality / 按质量过滤VCF文件"""
        qual_threshold = self.config.get('qual_threshold', 20)
        output_dir = self.config['output_dir']
        filtered_vcf = output_dir / f"{basename}.filtered.gz"
        
        self.logger.info(f"Filtering VCF file (QUAL >= {qual_threshold})...")
        self.logger.info(f"正在过滤VCF文件 (QUAL >= {qual_threshold})...")
        
        command = (
            f"bcftools filter -i 'QUAL>={qual_threshold}' "
            f"{vcf_file} -O z -o {filtered_vcf}"
        )
        
        try:
            self._execute_command(
                command,
                f"VCF quality filtering / VCF质量过滤 (QUAL>={qual_threshold})"
            )
            return filtered_vcf
        except ANNOVARCommandError:
            return None
    
    def _convert_vcf_to_annovar_format(self, vcf_file: Path, basename: str) -> bool:
        """Convert VCF to ANNOVAR format / 转换VCF为ANNOVAR格式"""
        output_dir = self.config['output_dir']
        annovar_path = self.config['annovar_path']
        annovar_vcf = output_dir / f"{basename}.annovar.vcf"
        
        self.logger.info("Converting VCF to ANNOVAR format...")
        self.logger.info("正在转换VCF为ANNOVAR格式...")
        
        command = (
            f"perl {annovar_path}/convert2annovar.pl "
            f"-format vcf4 {vcf_file} > {annovar_vcf}"
        )
        
        try:
            self._execute_command(
                command,
                "VCF to ANNOVAR format conversion / VCF到ANNOVAR格式转换"
            )
            self._intermediate_files['annovar_vcf'] = annovar_vcf
            self._intermediate_files['vcf_basename'] = basename
            return True
        except ANNOVARCommandError:
            return False
    
    def step4_annotate_variants(self) -> bool:
        """
        Step 4: Annotate variants using ANNOVAR.
        步骤4: 使用ANNOVAR注释变异。
        
        Returns:
            bool: Success status / 成功状态
        """
        if 'annovar_vcf' not in self._intermediate_files:
            raise ANNOVARError(
                "ANNOVAR VCF file not available. Please run step 3 first.\n"
                "ANNOVAR VCF文件不可用。请先运行步骤3。"
            )
        
        annovar_vcf = self._intermediate_files['annovar_vcf']
        annovar_path = self.config['annovar_path']
        database_path = self.config['database_path']
        build_ver = self.config['build_ver']
        output_dir = self.config['output_dir']
        vcf_basename = self._intermediate_files['vcf_basename']
        output_prefix = output_dir / vcf_basename
        
        self.logger.info("Annotating variants with ANNOVAR...")
        self.logger.info("正在使用ANNOVAR注释变异...")
        
        command = (
            f"perl {annovar_path}/annotate_variation.pl "
            f"--geneanno -dbtype refGene --buildver {build_ver} "
            f"{annovar_vcf} {database_path} -out {output_prefix}"
        )
        
        try:
            self._execute_command(
                command,
                "Variant annotation with ANNOVAR / 使用ANNOVAR进行变异注释"
            )
            
            # Collect output files / 收集输出文件
            self._collect_output_files(output_prefix)
            return True
            
        except ANNOVARCommandError:
            return False
    
    def _collect_output_files(self, output_prefix: Path) -> None:
        """Collect and log output files / 收集并记录输出文件"""
        expected_files = [
            output_prefix.with_suffix('.exonic_variant_function'),
            output_prefix.with_suffix('.variant_function'),
            output_prefix.with_suffix('.log')
        ]
        
        self._output_files = [f for f in expected_files if f.exists()]
        
        self.logger.info("Annotation completed successfully!")
        self.logger.info("注释成功完成!")
        self.logger.info("Output files / 输出文件:")
        
        for output_file in self._output_files:
            file_size = output_file.stat().st_size
            self.logger.info(f"  - {output_file.name} ({file_size} bytes)")
    
    def run_full_pipeline(self) -> bool:
        """
        Execute the complete ANNOVAR annotation pipeline.
        执行完整的ANNOVAR注释流水线。
        
        Returns:
            bool: Overall success status / 整体成功状态
        """
        self.logger.info("Starting ANNOVAR annotation pipeline...")
        self.logger.info("开始ANNOVAR注释流水线...")
        
        pipeline_steps = [
            (self.step1_convert_gff3_to_genepred, "GFF3 to GenePred conversion / GFF3到GenePred转换"),
            (self.step2_extract_transcript_sequences, "Transcript sequence extraction / 转录本序列提取"),
            (self.step3_process_and_convert_vcf, "VCF processing and conversion / VCF处理和转换"),
            (self.step4_annotate_variants, "Variant annotation / 变异注释")
        ]
        
        for step_num, (step_func, step_desc) in enumerate(pipeline_steps, 1):
            self.logger.info(f"\n=== Step {step_num}: {step_desc} ===")
            self.logger.info(f"=== 步骤 {step_num}: {step_desc} ===")
            
            try:
                if not step_func():
                    self.logger.error(f"Pipeline failed at step {step_num}: {step_desc}")
                    self.logger.error(f"流水线在步骤 {step_num} 失败: {step_desc}")
                    return False
                
                self.logger.info(f"Step {step_num} completed successfully.")
                self.logger.info(f"步骤 {step_num} 成功完成。")
                
            except Exception as e:
                self.logger.error(f"Unexpected error in step {step_num}: {e}")
                self.logger.error(f"步骤 {step_num} 发生意外错误: {e}")
                return False
        
        # Generate summary report / 生成总结报告
        self._generate_summary_report()
        
        self.logger.info("\n=== Pipeline completed successfully! ===")
        self.logger.info("=== 流水线成功完成! ===")
        
        return True
    
    def run_single_step(self, step_number: int) -> bool:
        """
        Execute a single step of the pipeline.
        执行流水线的单个步骤。
        
        Args:
            step_number (int): Step number (1-4) / 步骤编号（1-4）
            
        Returns:
            bool: Step success status / 步骤成功状态
            
        Raises:
            ANNOVARValidationError: If step number is invalid / 如果步骤编号无效
        """
        step_functions = {
            1: self.step1_convert_gff3_to_genepred,
            2: self.step2_extract_transcript_sequences,
            3: self.step3_process_and_convert_vcf,
            4: self.step4_annotate_variants
        }
        
        if step_number not in step_functions:
            raise ANNOVARValidationError(
                f"Invalid step number: {step_number}. Must be 1-4.\n"
                f"无效的步骤编号: {step_number}。必须是1-4。"
            )
        
        step_func = step_functions[step_number]
        step_name = step_func.__name__
        
        self.logger.info(f"Executing single step {step_number}: {step_name}")
        self.logger.info(f"执行单个步骤 {step_number}: {step_name}")
        
        return step_func()
    
    def get_output_files(self) -> List[Path]:
        """
        Get list of generated output files.
        获取生成的输出文件列表。
        
        Returns:
            List[Path]: List of output file paths / 输出文件路径列表
        """
        return self._output_files.copy()
    
    def get_intermediate_files(self) -> Dict[str, Path]:
        """
        Get dictionary of intermediate files.
        获取中间文件字典。
        
        Returns:
            Dict[str, Path]: Dictionary of intermediate file paths / 中间文件路径字典
        """
        return self._intermediate_files.copy()
    
    def _generate_summary_report(self) -> None:
        """Generate comprehensive summary report / 生成综合总结报告"""
        report_file = self.config['output_dir'] / "annovar_annotation_summary.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("ANNOVAR Annotation Summary Report\n")
            f.write("ANNOVAR注释总结报告\n")
            f.write("=" * 60 + "\n\n")
            
            # Input files section / 输入文件部分
            f.write("Input Files / 输入文件:\n")
            f.write("-" * 30 + "\n")
            f.write(f"GFF3 file / GFF3文件: {self.config['gff3_file']}\n")
            f.write(f"Genome file / 基因组文件: {self.config['genome_file']}\n")
            f.write(f"VCF file / VCF文件: {self.config['vcf_file']}\n\n")
            
            # Configuration section / 配置部分
            f.write("Configuration / 配置:\n")
            f.write("-" * 30 + "\n")
            f.write(f"Build version / 构建版本: {self.config['build_ver']}\n")
            f.write(f"ANNOVAR path / ANNOVAR路径: {self.config['annovar_path']}\n")
            f.write(f"Database path / 数据库路径: {self.config['database_path']}\n")
            f.write(f"Output directory / 输出目录: {self.config['output_dir']}\n")
            f.write(f"Quality threshold / 质量阈值: {self.config.get('qual_threshold', 20)}\n")
            f.write(f"Skip VCF filter / 跳过VCF过滤: {'Yes/是' if self.config.get('skip_vcf_filter', True) else 'No/否'}\n")
            f.write(f"Skip GFF fix / 跳过GFF修复: {'Yes/是' if self.config.get('skip_gff_fix', False) else 'No/否'}\n\n")
            
            # Output files section / 输出文件部分
            f.write("Output Files / 输出文件:\n")
            f.write("-" * 30 + "\n")
            if self._output_files:
                for i, output_file in enumerate(self._output_files, 1):
                    file_size = output_file.stat().st_size if output_file.exists() else 0
                    f.write(f"{i}. {output_file.name} ({file_size} bytes)\n")
                    f.write(f"   Full path / 完整路径: {output_file}\n")
            else:
                f.write("No output files generated / 未生成输出文件\n")
            
            # Intermediate files section / 中间文件部分
            if self._intermediate_files:
                f.write(f"\nIntermediate Files / 中间文件:\n")
                f.write("-" * 30 + "\n")
                for name, path in self._intermediate_files.items():
                    f.write(f"{name}: {path}\n")
        
        self.logger.info(f"Summary report generated: {report_file}")
        self.logger.info(f"总结报告已生成: {report_file}")
