# =============================================================================
# File: biohelpers/cli/annovar_cli.py
# =============================================================================
"""
ANNOVAR CLI Entry Point for biohelpers
biohelpers的ANNOVAR CLI入口点

Main CLI entry point for ANNOVAR annotation functionality.
ANNOVAR注释功能的主要CLI入口点。
"""

import sys
from ..annovar.cli import main as annovar_cli_main


def main():
    """
    Entry point for biohelpers ANNOVAR CLI.
    biohelpers ANNOVAR CLI的入口点。
    """
    return annovar_cli_main()


if __name__ == "__main__":
    sys.exit(main())