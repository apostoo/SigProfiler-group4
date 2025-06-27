#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:45:38 2019

@author: mishugeb
"""
from SigProfilerExtractor import sigpro as sig
import SigProfilerExtractor as spe_mod
import os


def run_matrix_96():
    data = sig.importdata("matrix")
    sig.sigProfilerExtractor(
        "matrix",
        "test_matrix_96_output",
        data,
        minimum_signatures=3,
        maximum_signatures=3,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
    )


def run_matrix_78():
    data = sig.importdata("matrix_DBS")
    sig.sigProfilerExtractor(
        "matrix",
        "test_matrix_78_output",
        data,
        minimum_signatures=3,
        maximum_signatures=3,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
    )


def run_matrix_83():
    data = sig.importdata("matrix_ID")
    sig.sigProfilerExtractor(
        "matrix",
        "test_matrix_83_output",
        data,
        exome=True,
        reference_genome="GRCh38",
        minimum_signatures=3,
        maximum_signatures=3,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
    )


def run_vcf():
    vcf_data = sig.importdata("vcf")
    sig.sigProfilerExtractor(
        "vcf",
        "test_vcf_output",
        vcf_data,
        minimum_signatures=3,
        maximum_signatures=3,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
        use_ssnmf=True,
        modelNum=1,
    )


def run_matrix_48():
    data = sig.importdata("matrix_CNV")
    sig.sigProfilerExtractor(
        "matrix",
        "test_matrix_48_output",
        data,
        minimum_signatures=3,
        maximum_signatures=3,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
    )


def run_seg_48():
    data = sig.importdata("seg:BATTENBERG")
    sig.sigProfilerExtractor(
        "seg:BATTENBERG",
        "test_segCNV_output",
        data,
        minimum_signatures=3,
        maximum_signatures=3,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
    )


def run_matrix_32():
    data = sig.importdata("matrix_SV")
    sig.sigProfilerExtractor(
        "matrix",
        "test_matrix_32_output",
        data,
        minimum_signatures=3,
        maximum_signatures=3,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
    )


def run_matobj():
    data = sig.importdata("matobj")
    sig.sigProfilerExtractor(
        "matobj",
        "test_matobj_output",
        data,
        minimum_signatures=3,
        maximum_signatures=3,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
    )


def run_csv():
    data = sig.importdata("csv")
    sig.sigProfilerExtractor(
        "csv",
        "test_csv_output",
        data,
        minimum_signatures=3,
        maximum_signatures=3,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
    )


def run_supervised_vcf():
    """
    Test the VCF processing functionality of supervised SSNMF.
    It will generate the Y matrix and save it to the Y_matrices directory.
    """
    vcf_data = sig.importdata("vcf")
    sig.sigProfilerExtractor(
        "vcf",
        "test_supervised_vcf_output",
        vcf_data,
        minimum_signatures=4,
        maximum_signatures=4,
        nmf_replicates=50,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
        # Supervised SSNMF parameters
        use_ssnmf=True,
        modelNum=3,  # Frobenius + Frobenius supervised model
        label_pattern="first_two_chars",  # Extract labels from file name prefixes
        lam=0.1,  # Supervised learning weight parameter, controlling the importance of the label term in the objective function
        # Optional: custom label mapping
        # custom_label_mapping={
        #     "PD3851a.vcf": "BRCA",
        #     "PD3890a.vcf": "LUAD",
        #     "PD3904a.vcf": "LUAD",
        #     "PD3905a.vcf": "BRCA",
        #     "PD3945a.vcf": "BRCA"
        # }
    )



if __name__ == "__main__":
    # run_matrix_96()
    # run_matrix_78()
    # run_matrix_83()
    # run_matrix_48()
    # run_matrix_32()
    # run_seg_48()
    # run_vcf()  # Modified VCF test (unsupervised SSNMF)
    run_supervised_vcf()  # Supervised SSNMF test (requires Y matrix)
    # run_matobj()
    # run_csv()
