# =============================================================================
# File: biohelpers/annovar/exceptions.py
# =============================================================================
"""
ANNOVAR-specific Exception Classes
ANNOVAR专用异常类

Custom exception classes for ANNOVAR annotation workflow.
ANNOVAR注释工作流程的自定义异常类。
"""

from typing import Optional


class ANNOVARError(Exception):
    """
    Base exception for ANNOVAR-related errors.
    ANNOVAR相关错误的基础异常类。
    
    Args:
        message (str): Error message / 错误消息
        error_code (Optional[int]): Error code / 错误代码
    """
    
    def __init__(self, message: str, error_code: Optional[int] = None):
        super().__init__(message)
        self.message = message
        self.error_code = error_code
    
    def __str__(self):
        if self.error_code:
            return f"[ANNOVAR Error {self.error_code}] {self.message}"
        return f"[ANNOVAR Error] {self.message}"


class ANNOVARValidationError(ANNOVARError):
    """
    Exception for validation errors in ANNOVAR workflow.
    ANNOVAR工作流程中验证错误的异常。
    """
    pass


class ANNOVARFileNotFoundError(ANNOVARError):
    """
    Exception for missing files in ANNOVAR workflow.
    ANNOVAR工作流程中文件缺失的异常。
    """
    
    def __init__(self, file_path: str, file_type: str = "file"):
        message = (
            f"{file_type.title()} not found: {file_path}\n"
            f"{file_type}未找到: {file_path}"
        )
        super().__init__(message)
        self.file_path = file_path
        self.file_type = file_type


class ANNOVARCommandError(ANNOVARError):
    """
    Exception for command execution failures in ANNOVAR workflow.
    ANNOVAR工作流程中命令执行失败的异常。
    """
    
    def __init__(self, message: str, command: Optional[str] = None, 
                 return_code: Optional[int] = None, stderr: Optional[str] = None):
        super().__init__(message)
        self.command = command
        self.return_code = return_code
        self.stderr = stderr
