# =============================================================================
# File: biohelpers/cli/__init__.py
# =============================================================================
"""
Command Line Interface Package for biohelpers
biohelpers的命令行接口包

Centralized CLI management for all biohelpers modules.
所有biohelpers模块的集中CLI管理。
"""

from .annovar_cli import main as annovar_main

__all__ = ['annovar_main']