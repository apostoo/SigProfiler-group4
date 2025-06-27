#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MC3 to VCF Converter
A tool to convert TCGA MC3 format files to standard VCF format.

Author: AI Assistant
Date: 2024
"""

import pandas as pd
import os
import numpy as np
from pathlib import Path
import argparse
import sys

class MC3ToVCFConverter:
    """
    MC3 to VCF Converter class.
    """
    
    def __init__(self, reference_genome="GRCh37"):
        self.reference_genome = reference_genome
        
    def convert_mc3_to_vcf(self, mc3_file, output_dir, project_name="MC3_Project"):
        """
        Converts an MC3 file to VCF format.
        
        Parameters:
        mc3_file (str): Path to the MC3 file.
        output_dir (str): Output directory.
        project_name (str): Project name.
        
        Returns:
        str: Path to the output directory.
        """
        print(f"=== Starting MC3 to VCF Conversion ===")
        print(f"Input file: {mc3_file}")
        print(f"Output directory: {output_dir}")
        print(f"Reference genome: {self.reference_genome}")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Read MC3 file
        try:
            mc3_data = pd.read_csv(mc3_file, sep='\t')
            print(f"Successfully read MC3 file: {len(mc3_data)} rows.")
        except Exception as e:
            print(f"Error: Failed to read MC3 file: {e}")
            return None
        
        # Validate required columns
        required_columns = ['sample', 'chr', 'start', 'reference', 'alt']
        missing_columns = [col for col in required_columns if col not in mc3_data.columns]
        if missing_columns:
            print(f"Error: Missing required columns: {missing_columns}")
            print(f"Available columns: {list(mc3_data.columns)}")
            return None
        
        # Clean data
        mc3_data = self._clean_mc3_data(mc3_data)
        
        # Group by sample
        samples = mc3_data['sample'].unique()
        print(f"Found {len(samples)} samples: {samples[:5]}{'...' if len(samples) > 5 else ''}")
        
        # Conversion statistics
        total_variants = 0
        converted_variants = 0
        
        # Process each sample
        for sample in samples:
            sample_data = mc3_data[mc3_data['sample'] == sample]
            vcf_file = os.path.join(output_dir, f"{sample}.vcf")
            
            print(f"\nProcessing sample: {sample} ({len(sample_data)} variants)")
            
            # Generate VCF file
            sample_converted = self._generate_vcf_for_sample(sample_data, sample, vcf_file)
            
            total_variants += len(sample_data)
            converted_variants += sample_converted
            
            print(f"Sample {sample}: Converted {sample_converted}/{len(sample_data)} variants.")
        
        # Generate summary
        self._generate_summary(output_dir, project_name, total_variants, converted_variants, samples)
        
        print(f"\n=== Conversion Complete ===")
        print(f"Total variants: {total_variants}")
        print(f"Successfully converted: {converted_variants}")
        if total_variants > 0:
            print(f"Conversion rate: {converted_variants/total_variants*100:.1f}%")
        print(f"Output written to: {output_dir}")
        
        return output_dir
    
    def _clean_mc3_data(self, data):
        """
        Cleans the MC3 data.
        
        Parameters:
        data (DataFrame): The MC3 DataFrame.
        
        Returns:
        DataFrame: The cleaned DataFrame.
        """
        print("Cleaning MC3 data...")
        
        # Remove empty rows
        original_count = len(data)
        data = data.dropna(subset=['sample', 'chr', 'start', 'reference', 'alt'])
        
        # Process chromosome names
        data['chr'] = data['chr'].astype(str).str.replace('chr', '')
        
        # Process position information
        data['start'] = pd.to_numeric(data['start'], errors='coerce')
        data = data.dropna(subset=['start'])
        
        # Process reference and alternate alleles
        data['reference'] = data['reference'].fillna('N')
        data['alt'] = data['alt'].fillna('N')
        
        # Remove invalid variants
        invalid_mask = (data['reference'] == 'N') & (data['alt'] == 'N')
        data = data[~invalid_mask]
        
        print(f"Data cleaning complete: {original_count} -> {len(data)} valid variants remaining.")
        return data
    
    def _generate_vcf_for_sample(self, sample_data, sample_name, vcf_file):
        """
        Generates a VCF file for a single sample.
        
        Parameters:
        sample_data (DataFrame): Data for the sample.
        sample_name (str): Name of the sample.
        vcf_file (str): Path to the output VCF file.
        
        Returns:
        int: Number of successfully converted variants.
        """
        converted_count = 0
        
        with open(vcf_file, 'w', encoding='utf-8') as f:
            # Write VCF header
            self._write_vcf_header(f)
            
            # Write data rows
            for _, row in sample_data.iterrows():
                try:
                    vcf_line = self._create_vcf_line(row)
                    if vcf_line:
                        f.write(vcf_line + '\n')
                        converted_count += 1
                except Exception as e:
                    print(f"Warning: Error processing variant: {e}")
                    continue
        
        return converted_count
    
    def _write_vcf_header(self, file_handle):
        """
        Writes the VCF header.
        
        Parameters:
        file_handle: File handle.
        """
        header_line = "#CHROM\tPOS\tFILTER\tREF\tALT"
        file_handle.write(header_line + '\n')
    
    def _create_vcf_line(self, row):
        """
        Creates a single VCF data line.
        
        Parameters:
        row (Series): Data row (Pandas Series).
        
        Returns:
        str: A VCF-formatted string, or None if invalid.
        """
        try:
            # Extract basic info
            chrom = str(row['chr'])
            pos = int(row['start'])
            ref = str(row['reference'])
            alt = str(row['alt'])
            
            # Validate data
            if pd.isna(ref) or pd.isna(alt) or ref == 'N' or alt == 'N':
                return None
            
            # Handle insertions/deletions
            if ref == '-' or ref == 'NA':
                ref = 'N'
            if alt == '-' or alt == 'NA':
                alt = 'N'
            
            # Skip invalid variants
            if ref == 'N' and alt == 'N':
                return None
            
            # Set filter status
            filter_status = "PASS"
            
            # Assemble VCF line
            vcf_line = f"{chrom}\t{pos}\t{filter_status}\t{ref}\t{alt}"
            
            return vcf_line
            
        except Exception as e:
            print(f"Warning: Error creating VCF line: {e}")
            return None
    
    def _generate_summary(self, output_dir, project_name, total_variants, converted_variants, samples):
        """
        Generates a conversion summary file.
        
        Parameters:
        output_dir (str): Output directory.
        project_name (str): Project name.
        total_variants (int): Total variants.
        converted_variants (int): Converted variants.
        samples (list): List of samples.
        """
        summary_file = os.path.join(output_dir, "conversion_summary.txt")
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("MC3 to VCF Conversion Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"Project Name: {project_name}\n")
            f.write(f"Reference Genome: {self.reference_genome}\n")
            f.write(f"Conversion Time: {pd.Timestamp.now()}\n")
            f.write(f"Number of Samples: {len(samples)}\n")
            f.write(f"Total Variants: {total_variants}\n")
            f.write(f"Successfully Converted: {converted_variants}\n")
            if total_variants > 0:
                f.write(f"Conversion Rate: {converted_variants/total_variants*100:.1f}%\n")
            f.write(f"Output Directory: {output_dir}\n")
            f.write("\nSample List:\n")
            for sample in samples:
                f.write(f"  - {sample}.vcf\n")
        
        print(f"Conversion summary saved to: {summary_file}")

def main():
    """
    Main function.
    """
    parser = argparse.ArgumentParser(description='MC3 to VCF Converter')
    parser.add_argument('input_file', help='Path to the MC3 file')
    parser.add_argument('output_dir', help='Output directory')
    parser.add_argument('--reference', default='GRCh37', help='Reference genome version')
    parser.add_argument('--project', default='MC3_Project', help='Project name')
    parser.add_argument('--demo', action='store_true', help='Run in demo mode')
    
    args = parser.parse_args()
    
    if args.demo:
        # Demo Mode
        print("=== MC3 to VCF Converter Demo ===")
        
        # Create demo MC3 file
        demo_mc3 = "demo_mc3.txt"
        demo_data = """sample	chr	start	end	reference	alt	gene	effect	Amino_Acid_Change	DNA_VAF	SIFT	PolyPhen
TCGA-3C-AAAU-01	10	122668955	122668955	G	A	WDR11	3'UTR		0.39		
TCGA-3C-AAAU-01	10	8115874	8115875	-	A	GATA3	Frame_Shift_Ins	p.P409Afs*99	0.34		
TCGA-3C-AAAU-01	11	65272906	65272908	AAA	-	MALAT1	RNA		0.27		
TCGA-3C-AAAU-01	11	66082467	66082467	C	T	CD248	Missense_Mutation	p.E678K	0.07	tolerated(0.12)	benign(0.001)
TCGA-3C-AAAU-01	11	66193652	66193652	G	C	NPAS4	3'UTR		0.20		
TCGA-3C-AAAU-01	11	89896900	89896900	G	A	NAALAD2	Intron		0.51		
TCGA-3C-AAAU-01	12	56709872	56709872	C	G	CNPY2	5'UTR		0.26		
TCGA-3C-AAAU-01	14	37838785	37838785	A	G	MIPOL1	Missense_Mutation	p.M298V	0.81	tolerated(0.41)	benign(0.002)
TCGA-3C-AAAU-01	15	72039337	72039337	C	G	THSD4	Missense_Mutation	p.L733V	0.31	deleterious(0.04)	benign(0.444)
TCGA-3C-AAAU-01	15	75130677	75130677	G	A	ULK3	Silent	p.A392A	0.21		"""
        
        with open(demo_mc3, 'w') as f:
            f.write(demo_data)
        
        print(f"Demo MC3 file created: {demo_mc3}")
        
        # Convert demo file
        converter = MC3ToVCFConverter()
        converter.convert_mc3_to_vcf(demo_mc3, "demo_vcf_output", "Demo_MC3_Project")
        
        print("\nDemo complete! Check the files in the demo_vcf_output directory.")
        
    else:
        # Normal Conversion Mode
        if not os.path.exists(args.input_file):
            print(f"Error: Input file not found: {args.input_file}")
            return
        
        converter = MC3ToVCFConverter(reference_genome=args.reference)
        converter.convert_mc3_to_vcf(args.input_file, args.output_dir, args.project)

if __name__ == "__main__":
    # If no command-line arguments are provided, run the demo.
    if len(sys.argv) == 1:
        # Create a dummy parser to trigger the help message for the demo
        parser = argparse.ArgumentParser(description='MC3 to VCF Converter - Demo Mode')
        parser.add_argument('--demo', action='store_true', help='Run in demo mode')
        args = parser.parse_args(['--demo'])
        main()
    else:
        main() 