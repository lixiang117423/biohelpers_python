#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Augustus Gene Prediction Complete Pipeline
Author: Automated Augustus Training and Prediction Pipeline
Date: 2025

This script implements a complete Augustus gene prediction pipeline including:
1. Create new species model
2. Prepare training data
3. Model training and optimization
4. Test set prediction
5. Result evaluation and Excel report generation
6. Format conversion
"""

import argparse
import os
import sys
import subprocess
import re
import logging
import pandas as pd
from pathlib import Path
from datetime import datetime


class AugustusTrainer:
    """Augustus training and prediction pipeline manager"""
    
    def __init__(self, config):
        """Initialize configuration"""
        self.config = config
        self.setup_logging()
        self.validate_inputs()
        
    def setup_logging(self):
        """Setup logging configuration"""
        log_file = os.path.join(self.config['output_dir'], 'augustus_pipeline.log')
        os.makedirs(self.config['output_dir'], exist_ok=True)
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file, encoding='utf-8'),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def validate_inputs(self):
        """Validate input files and paths"""
        self.logger.info("Validating input parameters...")
        
        # Check Augustus path
        if not os.path.exists(self.config['augustus_path']):
            raise FileNotFoundError(f"Augustus path does not exist: {self.config['augustus_path']}")
            
        # Check input files
        for file_key in ['genome_file', 'gff_file']:
            if not os.path.exists(self.config[file_key]):
                raise FileNotFoundError(f"Input file does not exist: {self.config[file_key]}")
                
        # Create output directory
        os.makedirs(self.config['output_dir'], exist_ok=True)
        
        self.logger.info("Input validation completed")
        
    def run_command(self, command, description=""):
        """Execute shell command"""
        self.logger.info(f"Executing: {description}")
        self.logger.debug(f"Command: {command}")
        
        try:
            result = subprocess.run(
                command, 
                shell=True, 
                check=True, 
                capture_output=True, 
                text=True,
                encoding='utf-8'
            )
            if result.stdout:
                self.logger.debug(f"Output: {result.stdout}")
            return result
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command execution failed: {e}")
            self.logger.error(f"Error output: {e.stderr}")
            raise
            
    def step1_create_species(self):
        """Step 1: Create new species model"""
        self.logger.info("=" * 50)
        self.logger.info("Step 1: Creating new species model")
        
        new_species_script = os.path.join(self.config['augustus_path'], 'new_species.pl')
        command = f"perl {new_species_script} --species={self.config['species_name']}"
        
        try:
            self.run_command(command, "Create new species model")
            self.logger.info(f"Successfully created species model: {self.config['species_name']}")
        except subprocess.CalledProcessError:
            self.logger.warning("Species model may already exist, continuing...")
            
    def step2_prepare_training_data(self):
        """Step 2: Prepare training data"""
        self.logger.info("=" * 50)
        self.logger.info("Step 2: Preparing training data")
        
        gff2gb_script = os.path.join(self.config['augustus_path'], 'gff2gbSmallDNA.pl')
        output_file = os.path.join(self.config['output_dir'], 'training_set.gb')
        
        command = (f"perl {gff2gb_script} "
                  f"{self.config['gff_file']} "
                  f"{self.config['genome_file']} "
                  f"{self.config['flank_length']} "
                  f"{output_file}")
        
        self.run_command(command, "Generate Augustus training data")
        
        self.config['training_file'] = output_file
        self.logger.info(f"Training data generated: {output_file}")
        
    def step3_split_dataset(self):
        """Step 3: Split dataset"""
        self.logger.info("=" * 50)
        self.logger.info("Step 3: Splitting training and test sets")
        
        # Count total genes
        with open(self.config['training_file'], 'r') as f:
            content = f.read()
            total_genes = content.count('LOCUS')
            
        self.logger.info(f"Detected total genes: {total_genes}")
        
        if total_genes < 100:
            raise ValueError("Total genes less than 100, insufficient for splitting and evaluation")
            
        # Calculate test set size (randomSplit.pl puts the specified number in .train and rest in .test)
        # To get 80% training and 20% testing, we need to specify the training count
        train_count = int(total_genes * self.config['train_ratio'])
        test_count = total_genes - train_count
        
        self.logger.info(f"Total genes: {total_genes}")
        self.logger.info(f"Training genes: {train_count} ({self.config['train_ratio']*100:.1f}%)")
        self.logger.info(f"Testing genes: {test_count} ({(1-self.config['train_ratio'])*100:.1f}%)")
        
        # Split dataset - randomSplit.pl takes the number for the .train file
        random_split_script = os.path.join(self.config['augustus_path'], 'randomSplit.pl')
        command = f"perl {random_split_script} {self.config['training_file']} {train_count}"
        
        self.run_command(command, "Split dataset")
        
        self.config['train_file'] = self.config['training_file'] + '.train'
        self.config['test_file'] = self.config['training_file'] + '.test'
        
        # Verify the split by checking file sizes
        if os.path.exists(self.config['train_file']) and os.path.exists(self.config['test_file']):
            train_size = os.path.getsize(self.config['train_file'])
            test_size = os.path.getsize(self.config['test_file'])
            
            # Count genes in each file to verify
            with open(self.config['train_file'], 'r') as f:
                actual_train_genes = f.read().count('LOCUS')
            with open(self.config['test_file'], 'r') as f:
                actual_test_genes = f.read().count('LOCUS')
                
            self.logger.info(f"Verification - Training file: {actual_train_genes} genes ({train_size/1024/1024:.1f} MB)")
            self.logger.info(f"Verification - Testing file: {actual_test_genes} genes ({test_size/1024/1024:.1f} MB)")
            
            # Check if split is correct (training should be larger for ratio > 0.5)
            if self.config['train_ratio'] > 0.5 and train_size < test_size:
                self.logger.warning("WARNING: Training file is smaller than test file!")
                self.logger.warning("This might indicate an issue with randomSplit.pl behavior")
                self.logger.warning(f"Expected {train_count} training genes, got {actual_train_genes}")
                
                # If the split seems reversed, swap the files
                if actual_train_genes < actual_test_genes and self.config['train_ratio'] > 0.5:
                    self.logger.info("Swapping train and test files to correct the split...")
                    temp_file = self.config['train_file'] + '.temp'
                    os.rename(self.config['train_file'], temp_file)
                    os.rename(self.config['test_file'], self.config['train_file'])
                    os.rename(temp_file, self.config['test_file'])
                    
                    self.logger.info("Files swapped successfully")
                    self.logger.info(f"Final - Training: {actual_test_genes} genes")
                    self.logger.info(f"Final - Testing: {actual_train_genes} genes")
        
        self.logger.info("Dataset splitting completed")
        
    def step4_train_model(self):
        """Step 4: Train model"""
        self.logger.info("=" * 50)
        self.logger.info("Step 4: Training model")
        
        # etraining
        etraining_bin = os.path.join(self.config['augustus_path'], 'etraining')
        command = f"{etraining_bin} --species={self.config['species_name']} {self.config['train_file']}"
        
        self.run_command(command, "etraining parameter training")
        self.logger.info("etraining completed")
        
        # optimize_augustus
        optimize_script = os.path.join(self.config['augustus_path'], 'optimize_augustus.pl')
        optimize_log = os.path.join(self.config['output_dir'], f'optimize_{self.config["species_name"]}.log')
        
        command = f"perl {optimize_script} --species={self.config['species_name']} {self.config['test_file']} > {optimize_log} 2>&1"
        
        self.run_command(command, "Model parameter optimization")
        self.logger.info("Model optimization completed")
        
    def step5_predict_test_set(self):
        """Step 5: Predict test set"""
        self.logger.info("=" * 50)
        self.logger.info("Step 5: Predicting test set")
        
        augustus_bin = os.path.join(self.config['augustus_path'], 'augustus')
        prediction_file = os.path.join(self.config['output_dir'], 'prediction_result.gff')
        
        command = f"{augustus_bin} --species={self.config['species_name']} {self.config['test_file']} > {prediction_file}"
        
        self.run_command(command, "Predict test set")
        
        self.config['prediction_file'] = prediction_file
        self.logger.info(f"Prediction results saved: {prediction_file}")
        
    def step6_parse_evaluation_results(self):
        """Step 6: Parse evaluation results"""
        self.logger.info("=" * 50)
        self.logger.info("Step 6: Parsing evaluation results")
        
        with open(self.config['prediction_file'], 'r') as f:
            content = f.read()
            
        # Extract evaluation data
        evaluation_data = self.extract_evaluation_metrics(content)
        
        # Generate Excel report
        self.generate_excel_report(evaluation_data)
        
        return evaluation_data
        
    def extract_evaluation_metrics(self, content):
        """Extract evaluation metrics"""
        evaluation = {}
        
        # Extract nucleotide level sensitivity and specificity
        nucleotide_pattern = r'nucleotide level\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|'
        nucleotide_match = re.search(nucleotide_pattern, content)
        if nucleotide_match:
            evaluation['nucleotide_sensitivity'] = float(nucleotide_match.group(1))
            evaluation['nucleotide_specificity'] = float(nucleotide_match.group(2))
            
        # Extract exon level data
        exon_pattern = r'exon level\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|.*?\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|'
        exon_match = re.search(exon_pattern, content, re.DOTALL)
        if exon_match:
            evaluation['exon_pred_total'] = int(exon_match.group(1))
            evaluation['exon_anno_total'] = int(exon_match.group(2))
            evaluation['exon_tp'] = int(exon_match.group(3))
            evaluation['exon_sensitivity'] = float(exon_match.group(4))
            evaluation['exon_specificity'] = float(exon_match.group(5))
            
        # Extract gene level data
        gene_pattern = r'gene level\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|'
        gene_match = re.search(gene_pattern, content)
        if gene_match:
            evaluation['gene_pred'] = int(gene_match.group(1))
            evaluation['gene_anno'] = int(gene_match.group(2))
            evaluation['gene_tp'] = int(gene_match.group(3))
            evaluation['gene_fp'] = int(gene_match.group(4))
            evaluation['gene_fn'] = int(gene_match.group(5))
            evaluation['gene_sensitivity'] = float(gene_match.group(6))
            evaluation['gene_specificity'] = float(gene_match.group(7))
            
        return evaluation
        
    def generate_excel_report(self, evaluation_data):
        """Generate Excel evaluation report with bilingual support"""
        self.logger.info("Generating Excel evaluation report")
        
        # Create evaluation results data
        results_data_en = []
        results_data_zh = []
        
        # Nucleotide level
        if 'nucleotide_sensitivity' in evaluation_data:
            results_data_en.extend([
                {
                    'Evaluation Level': 'Nucleotide Level',
                    'Metric': 'Sensitivity',
                    'Value': evaluation_data['nucleotide_sensitivity'],
                    'Description': 'Proportion of correctly predicted nucleotides, reflects model ability to find true genes'
                },
                {
                    'Evaluation Level': 'Nucleotide Level',
                    'Metric': 'Specificity',
                    'Value': evaluation_data['nucleotide_specificity'],
                    'Description': 'Proportion of accurately predicted nucleotides, reflects model prediction precision'
                }
            ])
            
            results_data_zh.extend([
                {
                    '评估级别': '核苷酸水平',
                    '评估指标': '敏感性',
                    '数值': evaluation_data['nucleotide_sensitivity'],
                    '说明': '正确预测的核苷酸比例，反映模型找到真实基因的能力'
                },
                {
                    '评估级别': '核苷酸水平',
                    '评估指标': '特异性',
                    '数值': evaluation_data['nucleotide_specificity'],
                    '说明': '预测准确的核苷酸比例，反映模型预测精度'
                }
            ])
            
        # Exon level
        if 'exon_sensitivity' in evaluation_data:
            results_data_en.extend([
                {
                    'Evaluation Level': 'Exon Level',
                    'Metric': 'Total Predicted Exons',
                    'Value': evaluation_data['exon_pred_total'],
                    'Description': 'Total number of exons predicted by the model'
                },
                {
                    'Evaluation Level': 'Exon Level',
                    'Metric': 'Total Annotated Exons',
                    'Value': evaluation_data['exon_anno_total'],
                    'Description': 'Total number of exons in reference annotation'
                },
                {
                    'Evaluation Level': 'Exon Level',
                    'Metric': 'True Positives',
                    'Value': evaluation_data['exon_tp'],
                    'Description': 'Number of correctly predicted exons (True Positive)'
                },
                {
                    'Evaluation Level': 'Exon Level',
                    'Metric': 'Sensitivity',
                    'Value': evaluation_data['exon_sensitivity'],
                    'Description': 'Proportion of correctly predicted exons among true exons'
                },
                {
                    'Evaluation Level': 'Exon Level',
                    'Metric': 'Specificity',
                    'Value': evaluation_data['exon_specificity'],
                    'Description': 'Proportion of correct exons among predicted exons'
                }
            ])
            
            results_data_zh.extend([
                {
                    '评估级别': '外显子水平',
                    '评估指标': '预测外显子总数',
                    '数值': evaluation_data['exon_pred_total'],
                    '说明': '模型预测的外显子总数量'
                },
                {
                    '评估级别': '外显子水平',
                    '评估指标': '注释外显子总数',
                    '数值': evaluation_data['exon_anno_total'],
                    '说明': '参考注释中的外显子总数量'
                },
                {
                    '评估级别': '外显子水平',
                    '评估指标': '正确预测数',
                    '数值': evaluation_data['exon_tp'],
                    '说明': '预测正确的外显子数量(True Positive)'
                },
                {
                    '评估级别': '外显子水平',
                    '评估指标': '敏感性',
                    '数值': evaluation_data['exon_sensitivity'],
                    '说明': '正确预测的外显子占真实外显子的比例'
                },
                {
                    '评估级别': '外显子水平',
                    '评估指标': '特异性',
                    '数值': evaluation_data['exon_specificity'],
                    '说明': '预测的外显子中正确的比例'
                }
            ])
            
        # Gene level
        if 'gene_sensitivity' in evaluation_data:
            results_data_en.extend([
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'Predicted Genes',
                    'Value': evaluation_data['gene_pred'],
                    'Description': 'Total number of genes predicted by the model'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'Annotated Genes',
                    'Value': evaluation_data['gene_anno'],
                    'Description': 'Total number of genes in reference annotation'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'True Positives (TP)',
                    'Value': evaluation_data['gene_tp'],
                    'Description': 'Number of completely correctly predicted genes'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'False Positives (FP)',
                    'Value': evaluation_data['gene_fp'],
                    'Description': 'Number of incorrectly predicted genes'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'False Negatives (FN)',
                    'Value': evaluation_data['gene_fn'],
                    'Description': 'Number of missed true genes'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'Sensitivity',
                    'Value': evaluation_data['gene_sensitivity'],
                    'Description': 'Proportion of correctly predicted genes among true genes'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'Specificity',
                    'Value': evaluation_data['gene_specificity'],
                    'Description': 'Proportion of correct genes among predicted genes'
                }
            ])
            
            results_data_zh.extend([
                {
                    '评估级别': '基因水平',
                    '评估指标': '预测基因数',
                    '数值': evaluation_data['gene_pred'],
                    '说明': '模型预测的基因总数'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '注释基因数',
                    '数值': evaluation_data['gene_anno'],
                    '说明': '参考注释中的基因总数'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '真阳性(TP)',
                    '数值': evaluation_data['gene_tp'],
                    '说明': '完全正确预测的基因数量'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '假阳性(FP)',
                    '数值': evaluation_data['gene_fp'],
                    '说明': '错误预测的基因数量'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '假阴性(FN)',
                    '数值': evaluation_data['gene_fn'],
                    '说明': '漏掉的真实基因数量'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '敏感性',
                    '数值': evaluation_data['gene_sensitivity'],
                    '说明': '正确预测的基因占真实基因的比例'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '特异性',
                    '数值': evaluation_data['gene_specificity'],
                    '说明': '预测基因中正确的比例'
                }
            ])
            
        # Create DataFrames
        df_results_en = pd.DataFrame(results_data_en)
        df_results_zh = pd.DataFrame(results_data_zh)
        
        # Create configuration DataFrames
        config_data_en = [
            ['Species Name', self.config['species_name']],
            ['Genome File', self.config['genome_file']],
            ['Annotation File', self.config['gff_file']],
            ['Training Ratio', f"{self.config['train_ratio']*100}%"],
            ['Flank Length', f"{self.config['flank_length']} bp"],
            ['Output Directory', self.config['output_dir']],
            ['Generation Time', datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
        ]
        df_config_en = pd.DataFrame(config_data_en, columns=['Parameter', 'Value'])
        
        config_data_zh = [
            ['物种名称', self.config['species_name']],
            ['基因组文件', self.config['genome_file']],
            ['注释文件', self.config['gff_file']],
            ['训练集比例', f"{self.config['train_ratio']*100}%"],
            ['侧翼长度', f"{self.config['flank_length']} bp"],
            ['输出目录', self.config['output_dir']],
            ['生成时间', datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
        ]
        df_config_zh = pd.DataFrame(config_data_zh, columns=['参数', '值'])
        
        # Save Excel files (English and Chinese versions)
        excel_file_en = os.path.join(self.config['output_dir'], 'augustus_evaluation_report_EN.xlsx')
        excel_file_zh = os.path.join(self.config['output_dir'], 'augustus_evaluation_report_ZH.xlsx')
        
        # English version
        with pd.ExcelWriter(excel_file_en, engine='openpyxl') as writer:
            df_config_en.to_excel(writer, sheet_name='Configuration', index=False)
            df_results_en.to_excel(writer, sheet_name='Evaluation Results', index=False)
            
            # Add explanation sheet
            explanation_data_en = [
                ['Term', 'Explanation'],
                ['Sensitivity', 'Also called recall, represents model ability to correctly identify true genes. Formula: TP/(TP+FN)'],
                ['Specificity', 'Represents model prediction accuracy. Formula: TP/(TP+FP)'],
                ['TP (True Positive)', 'Number of correctly predicted genes'],
                ['FP (False Positive)', 'Number of incorrectly predicted genes'],
                ['FN (False Negative)', 'Number of missed true genes'],
                ['Nucleotide Level', 'Prediction accuracy at DNA sequence base level'],
                ['Exon Level', 'Prediction accuracy at exon structure level'],
                ['Gene Level', 'Prediction accuracy at complete gene level'],
                ['Evaluation Suggestion', 'Generally, models with sensitivity>0.8 and specificity>0.8 are considered excellent']
            ]
            df_explanation_en = pd.DataFrame(explanation_data_en[1:], columns=explanation_data_en[0])
            df_explanation_en.to_excel(writer, sheet_name='Term Explanations', index=False)
            
        # Chinese version
        with pd.ExcelWriter(excel_file_zh, engine='openpyxl') as writer:
            df_config_zh.to_excel(writer, sheet_name='配置信息', index=False)
            df_results_zh.to_excel(writer, sheet_name='评估结果', index=False)
            
            # Add explanation sheet
            explanation_data_zh = [
                ['术语', '解释'],
                ['敏感性(Sensitivity)', '也称召回率，表示模型正确识别真实基因的能力，计算公式: TP/(TP+FN)'],
                ['特异性(Specificity)', '表示模型预测准确度，计算公式: TP/(TP+FP)'],
                ['TP (True Positive)', '真阳性，正确预测的基因数量'],
                ['FP (False Positive)', '假阳性，错误预测的基因数量'],
                ['FN (False Negative)', '假阴性，漏掉的真实基因数量'],
                ['核苷酸水平', '在DNA序列碱基层面的预测准确性'],
                ['外显子水平', '在外显子结构层面的预测准确性'],
                ['基因水平', '在完整基因层面的预测准确性'],
                ['评估建议', '一般认为敏感性>0.8、特异性>0.8的模型较为优秀']
            ]
            df_explanation_zh = pd.DataFrame(explanation_data_zh[1:], columns=explanation_data_zh[0])
            df_explanation_zh.to_excel(writer, sheet_name='术语解释', index=False)
            
        self.logger.info(f"Excel evaluation reports generated:")
        self.logger.info(f"  English version: {excel_file_en}")
        self.logger.info(f"  Chinese version: {excel_file_zh}")
        
    def step7_convert_to_gff3(self):
        """Step 7: Convert to GFF3 format"""
        self.logger.info("=" * 50)
        self.logger.info("Step 7: Converting to GFF3 format")
        
        gff3_file = os.path.join(self.config['output_dir'], 'prediction_result.gff3')
        
        # Use gffread to convert format
        command = f"gffread {self.config['prediction_file']} -o {gff3_file}"
        
        try:
            self.run_command(command, "Convert to GFF3 format")
            self.logger.info(f"GFF3 file generated: {gff3_file}")
        except subprocess.CalledProcessError:
            self.logger.warning("gffread conversion failed, attempting simple format conversion...")
            self.simple_gff_to_gff3_conversion(gff3_file)
            
    def simple_gff_to_gff3_conversion(self, output_file):
        """Simple GFF to GFF3 conversion"""
        with open(self.config['prediction_file'], 'r') as infile, \
             open(output_file, 'w') as outfile:
            
            outfile.write("##gff-version 3\n")
            
            for line in infile:
                if line.startswith('#') or line.strip() == '':
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    # Simple processing of attribute field to ensure GFF3 format compliance
                    attributes = fields[8]
                    if 'transcript_id' in attributes and 'gene_id' in attributes:
                        outfile.write(line)
                        
        self.logger.info(f"Simple format conversion completed: {output_file}")
        
    def run_complete_pipeline(self):
        """Run complete pipeline"""
        try:
            self.logger.info("Starting Augustus complete training and prediction pipeline")
            self.logger.info(f"Species name: {self.config['species_name']}")
            
            # Execute all steps
            self.step1_create_species()
            self.step2_prepare_training_data()
            self.step3_split_dataset()
            self.step4_train_model()
            self.step5_predict_test_set()
            evaluation_data = self.step6_parse_evaluation_results()
            self.step7_convert_to_gff3()
            
            self.logger.info("=" * 50)
            self.logger.info("🎉 Augustus pipeline execution completed!")
            self.logger.info(f"Result files saved in: {self.config['output_dir']}")
            
            # Print key evaluation results
            if evaluation_data:
                self.logger.info("\nKey evaluation results:")
                if 'nucleotide_sensitivity' in evaluation_data:
                    self.logger.info(f"  Nucleotide sensitivity: {evaluation_data['nucleotide_sensitivity']:.3f}")
                    self.logger.info(f"  Nucleotide specificity: {evaluation_data['nucleotide_specificity']:.3f}")
                if 'gene_sensitivity' in evaluation_data:
                    self.logger.info(f"  Gene sensitivity: {evaluation_data['gene_sensitivity']:.3f}")
                    self.logger.info(f"  Gene specificity: {evaluation_data['gene_specificity']:.3f}")
                    
        except Exception as e:
            self.logger.error(f"Pipeline execution failed: {e}")
            raise


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Augustus Gene Prediction Complete Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage Examples:
  python augustus_pipeline.py \\
    --species_name Rice_NLR_Model \\
    --genome_file genome.fa \\
    --gff_file annotations.gff3 \\
    --output_dir ./augustus_results \\
    --train_ratio 0.8 \\
    --flank_length 1000 \\
    --augustus_path /path/to/augustus/bin

Detailed Description:
  This script automatically executes the complete Augustus training and 
  prediction pipeline, including model training, parameter optimization, 
  prediction evaluation, and result report generation.
        """
    )
    
    # Required parameters
    parser.add_argument('--species_name', required=True,
                       help='New species model name (e.g., Rice_NLR_Model)')
    parser.add_argument('--genome_file', required=True,
                       help='Genome FASTA file path')
    parser.add_argument('--gff_file', required=True,
                       help='Gene annotation GFF3 file path')
    
    # Optional parameters
    parser.add_argument('--output_dir', default='./augustus_output',
                       help='Output directory path (default: ./augustus_output)')
    parser.add_argument('--augustus_path', 
                       default='/share/org/YZWL/yzwl_lixg/miniforge3/envs/Augustus_v.3.5.0/bin',
                       help='Augustus installation path')
    parser.add_argument('--train_ratio', type=float, default=0.8,
                       help='Training set ratio (default: 0.8)')
    parser.add_argument('--flank_length', type=int, default=1000,
                       help='Gene flanking length (default: 1000)')
    
    args = parser.parse_args()
    
    # Build configuration dictionary
    config = {
        'species_name': args.species_name,
        'genome_file': os.path.abspath(args.genome_file),
        'gff_file': os.path.abspath(args.gff_file),
        'output_dir': os.path.abspath(args.output_dir),
        'augustus_path': args.augustus_path,
        'train_ratio': args.train_ratio,
        'flank_length': args.flank_length
    }
    
    # Execute pipeline
    trainer = AugustusTrainer(config)
    trainer.run_complete_pipeline()


if __name__ == '__main__':
    main()