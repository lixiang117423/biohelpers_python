#!/usr/bin/env python3
"""
GBFF to GFF converter
从NCBI的GBFF文件中提取GFF格式的注释信息
"""

import re
import sys
from typing import Dict, List, Optional, Tuple


class GBFFParser:
    def __init__(self):
        self.gff_version = "##gff-version 3"
        self.sequence_region = None

    def parse_location(self, location_str: str) -> Tuple[int, int, str]:
        """
        解析GenBank位置字符串
        返回: (start, end, strand)
        """
        # 处理互补链
        strand = "+"
        if location_str.startswith("complement("):
            strand = "-"
            location_str = location_str[11:-1]  # 移除complement()

        # 处理join位置 - 暂时取第一个片段
        if location_str.startswith("join("):
            # 提取第一个位置
            match = re.search(r"(\d+)\.\.(\d+)", location_str)
            if match:
                start, end = int(match.group(1)), int(match.group(2))
            else:
                return 1, 1, strand
        else:
            # 简单位置 start..end
            if ".." in location_str:
                try:
                    start_str, end_str = location_str.split("..")
                    # 处理 <start 或 >end 的情况
                    start = int(re.sub(r"[<>]", "", start_str))
                    end = int(re.sub(r"[<>]", "", end_str))
                except ValueError:
                    return 1, 1, strand
            else:
                # 单个位置
                try:
                    pos = int(re.sub(r"[<>]", "", location_str))
                    start, end = pos, pos
                except ValueError:
                    return 1, 1, strand

        return start, end, strand

    def extract_attributes(self, feature_lines: List[str]) -> Dict[str, str]:
        """从feature行中提取属性"""
        attributes = {}
        current_qualifier = None
        current_value = ""

        for line in feature_lines:
            line = line.strip()
            if line.startswith("/"):
                # 保存上一个qualifier
                if current_qualifier:
                    attributes[current_qualifier] = current_value.strip('"')

                # 解析新的qualifier
                if "=" in line:
                    qualifier, value = line[1:].split("=", 1)
                    current_qualifier = qualifier
                    current_value = value
                else:
                    current_qualifier = line[1:]
                    current_value = ""
            else:
                # 继续上一个qualifier的值
                current_value += " " + line

        # 保存最后一个qualifier
        if current_qualifier:
            attributes[current_qualifier] = current_value.strip('"')

        return attributes

    def format_gff_attributes(
        self, attributes: Dict[str, str], feature_type: str
    ) -> str:
        """格式化GFF属性字段"""
        gff_attrs = []

        # ID和Name
        if "locus_tag" in attributes:
            gff_attrs.append(f"ID={attributes['locus_tag']}")
            gff_attrs.append(f"Name={attributes['locus_tag']}")
        elif "gene" in attributes:
            gff_attrs.append(f"ID={attributes['gene']}")
            gff_attrs.append(f"Name={attributes['gene']}")

        # 产品信息
        if "product" in attributes:
            # 清理产品名称中的特殊字符
            product = (
                attributes["product"]
                .replace(";", "%3B")
                .replace("=", "%3D")
                .replace("&", "%26")
            )
            gff_attrs.append(f"product={product}")

        # 基因信息
        if "gene" in attributes and feature_type != "gene":
            gff_attrs.append(f"gene={attributes['gene']}")

        # 蛋白质ID
        if "protein_id" in attributes:
            gff_attrs.append(f"protein_id={attributes['protein_id']}")

        # 翻译
        if "translation" in attributes:
            # 翻译序列通常很长，可以选择性包含
            translation = (
                attributes["translation"][:50] + "..."
                if len(attributes["translation"]) > 50
                else attributes["translation"]
            )
            gff_attrs.append(f"translation={translation}")

        # 注释
        if "note" in attributes:
            note = (
                attributes["note"]
                .replace(";", "%3B")
                .replace("=", "%3D")
                .replace("&", "%26")
            )
            gff_attrs.append(f"note={note}")

        return ";".join(gff_attrs) if gff_attrs else "."

    def parse_gbff_file(self, gbff_file: str, output_file: str):
        """解析GBFF文件并输出GFF格式"""
        with open(gbff_file, "r", encoding="utf-8") as f:
            content = f.read()

        with open(output_file, "w", encoding="utf-8") as out:
            # 写入GFF头部
            out.write(f"{self.gff_version}\n")

            # 按记录分割（每个//结束一个记录）
            records = content.split("//\n")

            for record in records:
                if not record.strip():
                    continue

                lines = record.split("\n")
                seqid = None
                seq_length = None
                in_features = False

                # 解析序列信息
                for line in lines:
                    if line.startswith("LOCUS"):
                        parts = line.split()
                        if len(parts) >= 3:
                            seqid = parts[1]
                            seq_length = parts[2]
                    elif line.startswith("ACCESSION"):
                        if not seqid:  # 如果LOCUS中没有合适的ID，使用ACCESSION
                            seqid = line.split()[1]

                if not seqid:
                    seqid = "unknown"

                # 写入序列区域信息
                if seq_length and seq_length.isdigit():
                    out.write(f"##sequence-region {seqid} 1 {seq_length}\n")

                # 解析features
                i = 0
                while i < len(lines):
                    line = lines[i]

                    if line.startswith("FEATURES"):
                        in_features = True
                        i += 1
                        continue

                    if line.startswith("ORIGIN") or line.startswith("//"):
                        in_features = False
                        break

                    if in_features and len(line) > 5 and line[5] != " ":
                        # 这是一个新的feature行
                        feature_match = re.match(r"     (\w+)\s+(.+)", line)
                        if feature_match:
                            feature_type = feature_match.group(1)
                            location = feature_match.group(2)

                            # 收集所有属于这个feature的行
                            feature_lines = []
                            i += 1
                            while i < len(lines) and (
                                lines[i].startswith("                     ")
                                or (len(lines[i]) > 20 and lines[i][:21].strip() == "")
                            ):
                                feature_lines.append(lines[i])
                                i += 1
                            i -= 1  # 回退一行

                            # 解析位置
                            start, end, strand = self.parse_location(location)

                            # 提取属性
                            attributes = self.extract_attributes(feature_lines)

                            # 格式化GFF行
                            source = "GenBank"
                            score = "."
                            phase = "." if feature_type != "CDS" else "0"

                            gff_attributes = self.format_gff_attributes(
                                attributes, feature_type
                            )

                            # 写入GFF行
                            gff_line = f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{gff_attributes}\n"
                            out.write(gff_line)

                    i += 1


def main():
    if len(sys.argv) != 3:
        print("使用方法: python gbff_to_gff.py <input.gbff> <output.gff>")
        print("示例: python gbff_to_gff.py sequence.gbff output.gff")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    parser = GBFFParser()

    try:
        parser.parse_gbff_file(input_file, output_file)
        print(f"成功转换: {input_file} -> {output_file}")
    except FileNotFoundError:
        print(f"错误: 找不到输入文件 {input_file}")
    except Exception as e:
        print(f"转换过程中出现错误: {e}")


if __name__ == "__main__":
    main()
