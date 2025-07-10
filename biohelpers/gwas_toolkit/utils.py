"""
GWAS分析工具函数模块 | GWAS Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional

class GWASLogger:
    """GWAS分析日志管理器 | GWAS Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "gwas_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir
    
    def run(self, cmd: str, description: str = "") -> str:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"开始 | Starting: {description}")
        
        cleaned_cmd = " ".join(cmd.strip().split())
        self.logger.info(f"执行命令 | Executing command: {cleaned_cmd}")
        
        try:
            result = subprocess.run(
                cleaned_cmd, 
                shell=True, 
                check=True, 
                capture_output=True, 
                text=True, 
                cwd=self.working_dir
            )
            
            if result.stdout.strip():
                self.logger.info(f"命令输出 | Command output:\n{result.stdout.strip()}")
            
            if result.stderr.strip():
                self.logger.warning(f"命令标准错误流(可能包含日志信息) | Command stderr (may contain log info):\n{result.stderr.strip()}")
            
            self.logger.info(f"完成 | Completed: {description}")
            return result.stdout
            
        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"命令执行失败 | Command execution failed: {cleaned_cmd}\n"
                f"返回码 | Return code: {e.returncode}\n"
                f"标准输出 | Stdout: {e.stdout.strip()}\n"
                f"标准错误 | Stderr: {e.stderr.strip()}"
            )
            raise

def check_dependencies():
    """检查依赖 | Check dependencies"""
    try:
        from statsmodels.stats.multitest import multipletests
    except ImportError:
        print("错误 | Error: 本脚本的FDR校正功能需要 `statsmodels` 库。| FDR correction requires `statsmodels` library.")
        print("请运行 | Please run: pip install statsmodels")
        sys.exit(1)
