from SigProfilerExtractor import sigpro as sig

from pathlib import Path

"""
    Generate exposure matrices.
    Generate one matrix per different signature numbers.
    For the project, I ran signatures from 2 to 15 due to time constraints.
"""
def generate_matrices(input_datatype, input_path):
    for i in range(7, 16):
        signatures_filename = f"{i}-signatures"
        
        output_file = f"../results/exposure-matrices/results-{signatures_filename}"

        sig.sigProfilerExtractor(input_datatype, output_file, input_path, minimum_signatures=i, maximum_signatures=i, nmf_replicates=30)


if __name__ == '__main__':
    generate_matrices("vcf", "../data/vcf_per_sample/")