#!/usr/bin/env python3
"""
FASTA文件字符清理脚本
删除FASTA序列数据中的特定字符，但保持序列名称行不变
"""

import argparse
import sys
import os


def clean_fasta(input_file, output_file, chars_to_remove):
    """
    清理FASTA文件中的特定字符
    
    参数:
    input_file: 输入FASTA文件路径
    output_file: 输出FASTA文件路径
    chars_to_remove: 要删除的字符字符串
    """
    try:
        with open(input_file, 'r', encoding='utf-8') as infile, \
             open(output_file, 'w', encoding='utf-8') as outfile:
            
            for line in infile:
                line = line.rstrip('\n\r')  # 移除行尾换行符
                
                # 如果是序列名称行（以>开头），直接写入，不做任何修改
                if line.startswith('>'):
                    outfile.write(line + '\n')
                else:
                    # 序列数据行，删除指定字符
                    cleaned_line = line
                    for char in chars_to_remove:
                        cleaned_line = cleaned_line.replace(char, '')
                    outfile.write(cleaned_line + '\n')
                    
        print(f"成功处理文件: {input_file} -> {output_file}")
        print(f"删除的字符: {chars_to_remove}")
        
    except FileNotFoundError:
        print(f"错误: 找不到输入文件 '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except PermissionError:
        print(f"错误: 没有权限访问文件", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"处理文件时发生错误: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(
        description='删除FASTA文件序列数据中的特定字符，保持序列名称不变',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python fasta_cleaner.py -i input.fasta -s ".*-" -o output.fasta
  python fasta_cleaner.py --input sequences.fa --string "N-" --output clean_sequences.fa
        """
    )
    
    # 添加命令行参数
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='输入FASTA文件路径'
    )
    
    parser.add_argument(
        '-s', '--string',
        required=True,
        help='要删除的字符串（如: ".*-"）'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='输出FASTA文件路径'
    )
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"错误: 输入文件 '{args.input}' 不存在", file=sys.stderr)
        sys.exit(1)
    
    # 检查输出目录是否存在
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        print(f"错误: 输出目录 '{output_dir}' 不存在", file=sys.stderr)
        sys.exit(1)
    
    # 执行清理操作
    clean_fasta(args.input, args.output, args.string)


if __name__ == '__main__':
    main()