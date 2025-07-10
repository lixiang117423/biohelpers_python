# biohelpers

生物信息学常用脚本工具包，旨在为生信分析提供便捷的Python函数与脚本集合。

---

## 快速开始

### 安装

```bash
# 推荐使用pip本地安装
pip install .
```

---

## 主要模块与用法

### 1. annovar 模块

#### 1.1 Python API 用法

```python
from biohelpers.annovar import ANNOVARAnnotator

config = {
    'gff3_file': 'genome.gff3',
    'genome_file': 'genome.fasta',
    'vcf_file': 'variants.vcf',
    'build_ver': 'OV',
    'annovar_path': '/path/to/annovar',  # ANNOVAR软件目录
    'database_path': './database',       # ANNOVAR数据库目录
    'output_dir': './annovar_output'     # 输出目录
}

annotator = ANNOVARAnnotator(config)
success = annotator.run_full_pipeline()
if success:
    print("注释完成！")
    print(annotator.get_output_files())
else:
    print("注释失败！")
```

#### 1.2 命令行用法

```bash
biohelpers-annovar \
  --gff3 genome.gff3 \
  --genome genome.fa \
  --vcf variants.vcf \
  --build-ver OV \
  --annovar-path /path/to/annovar \
  --database-path ./database \
  --output-dir ./annovar_output
```

**常用参数说明：**
- `--gff3`：GFF3注释文件路径
- `--genome`：基因组FASTA文件路径
- `--vcf`：VCF变异文件路径
- `--build-ver`：基因组构建版本（如hg38、OV等）
- `--annovar-path`：ANNOVAR软件目录
- `--database-path`：ANNOVAR数据库目录
- `--output-dir`：输出目录
- `--step`：仅执行指定步骤（1-4）

更多参数和示例请使用 `biohelpers-annovar --help` 查看。

---

### 2. cli 模块

cli模块用于统一管理命令行接口，当前主要入口为annovar相关命令。

---

### 3. comm 模块

comm模块为预留的通用功能扩展，目前暂无具体实现。

---

如需扩展新功能，请在对应目录下添加模块，并在此文档中补充用法说明。 