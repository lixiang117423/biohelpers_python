#!/usr/bin/env python3
"""
模块提取脚本 | Module Extraction Script
从打包文件中提取各个模块到对应目录 | Extract modules from packed file to corresponding directories

用法 | Usage: python extract_modules.py [打包文件路径 | packed_file_path] [目标目录 | target_directory]

Examples:
    python extract_modules.py gwas_modules.txt .
    python extract_modules.py tmp.txt /path/to/project
"""

import sys
import os
import re
from pathlib import Path

def extract_modules(packed_file_path, target_dir="."):
    """从打包文件中提取模块"""
    
    # 读取打包文件
    try:
        with open(packed_file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"错误 | Error: 找不到打包文件 | Cannot find packed file '{packed_file_path}'")
        sys.exit(1)
    except Exception as e:
        print(f"错误 | Error: 读取打包文件失败 | Failed to read packed file: {e}")
        sys.exit(1)
    
    # 使用正则表达式找到所有文件块
    file_pattern = r'# ===== FILE: (.+?) =====\n(.*?)(?=# ===== END FILE =====)'
    matches = re.findall(file_pattern, content, re.DOTALL)
    
    if not matches:
        print("错误 | Error: 在打包文件中未找到任何模块 | No modules found in packed file")
        sys.exit(1)
    
    target_path = Path(target_dir)
    extracted_files = []
    
    for file_path, file_content in matches:
        # 创建完整的文件路径
        full_path = target_path / file_path
        
        # 创建目录（如果不存在）
        full_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 写入文件
        try:
            with open(full_path, 'w', encoding='utf-8') as f:
                f.write(file_content.rstrip() + '\n')
            
            # 如果是Python脚本，设置执行权限
            if file_path.endswith('.py') and ('#!/usr/bin/env python3' in file_content or file_path == 'run_gwas.py'):
                os.chmod(full_path, 0o755)
            
            extracted_files.append(str(full_path))
            print(f"已提取 | Extracted: {full_path}")
            
        except Exception as e:
            print(f"错误 | Error: 写入文件失败 | Failed to write file '{full_path}': {e}")
    
    # 显示总结
    print(f"\n提取完成 | Extraction completed! 共提取了 | Total extracted {len(extracted_files)} 个文件 | files:")
    for file in extracted_files:
        print(f"  - {file}")
    
    # # 显示使用说明
    # print(f"\n使用说明 | Usage Instructions:")
    # print(f"1. 进入目标目录 | Enter target directory: cd {target_dir}")
    # print(f"2. 运行GWAS分析 | Run GWAS analysis:")
    # print(f"   python run_gwas.py --vcf your_file.vcf --pheno your_pheno.txt")
    # print(f"   或者 | or: python -m gwas_toolkit.main --vcf your_file.vcf --pheno your_pheno.txt")
    # print(f"\n项目结构 | Project Structure:")
    # print(f"  {target_dir}/")
    # print(f"  ├── gwas_toolkit/           # GWAS分析包 | GWAS analysis package")
    # print(f"  │   ├── __init__.py         # 包初始化 | Package initialization")
    # print(f"  │   ├── config.py           # 配置管理 | Configuration management")
    # print(f"  │   ├── utils.py            # 工具函数 | Utility functions")
    # print(f"  │   ├── data_processing.py  # 数据处理 | Data processing")
    # print(f"  │   ├── gwas_analysis.py    # GWAS分析核心 | GWAS analysis core")
    # print(f"  │   ├── results.py          # 结果处理 | Results processing")
    # print(f"  │   └── main.py             # 主程序 | Main program")
    # print(f"  └── run_gwas.py             # 简化的运行脚本 | Simplified runner script")

def main():
    """主函数 | Main function"""
    if len(sys.argv) < 2:
        print("用法 | Usage: python extract_modules.py <打包文件路径 | packed_file_path> [目标目录 | target_directory]")
        print("例如 | Examples:")
        print("  python extract_modules.py tmp.txt .")
        print("  python extract_modules.py gwas_modules.txt /path/to/project")
        sys.exit(1)
    
    packed_file = sys.argv[1]
    target_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    
    print(f"开始提取模块 | Starting module extraction...")
    print(f"打包文件 | Packed file: {packed_file}")
    print(f"目标目录 | Target directory: {target_dir}")
    print("-" * 50)
    
    extract_modules(packed_file, target_dir)

if __name__ == "__main__":
    main()