# biohelpers package
# __version__ = '0.1.0'

# Import core modules
from .parse_longest_mrna import *

# package/__init__.py
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("your-package")
except PackageNotFoundError:
    __version__ = "0.0.0.dev"  # 开发环境备用值