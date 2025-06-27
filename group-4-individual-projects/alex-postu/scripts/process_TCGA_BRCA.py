import pandas as pd 
from pathlib import Path

valid_bases = {'A', 'C', 'G', 'T'}
expected_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y']

# Check if a ref or alt base is valid.
def is_valid_base(base):
    return isinstance(base, str) and base in valid_bases

# Fix invalid bases. 
# Only works for cases when only one of the two is invalid.
def fix_invalid(row):
    ref, alt = row['reference'], row['alt']
    ref_valid = is_valid_base(ref)
    alt_valid = is_valid_base(alt)

    if not ref_valid and alt_valid:
        row['reference'] = alt
    elif ref_valid and not alt_valid:
        row['alt'] = ref
    return row


# Sort chromosomes in the right order
def chrom_order(chrom):
    if chrom.isdigit():
        return int(chrom)
    elif chrom == 'X':
        return 23
    elif chrom == 'Y':
        return 24
    else:
        return 25


def process_dataset(file):
    df = pd.read_csv(file, sep='\t')
    print(df)

    # Construct INFO field by joining relevant columns
    def build_info(row):
        info_fields = []
        for key in ['gene', 'effect', 'Amino_Acid_Change', 'DNA_VAF', 'SIFT', 'PolyPhen']:
            val = row[key]
            if pd.notna(val) and val != '':
                info_fields.append(f"{key}={val}")
        return ";".join(info_fields) if info_fields else '.'

    df['INFO'] = df.apply(build_info, axis=1)

    # Sort chromosomes
    chrom_order = {str(i): i for i in range(1, 23)}
    chrom_order.update({'X': 23, 'Y': 24, 'MT': 25})

    # Normalize chromosome column to string
    df['chr'] = df['chr'].astype(str)
    df['chr_sort'] = df['chr'].map(chrom_order)

    # Sort the dataframe before grouping
    df = df.sort_values(by=['sample', 'chr_sort', 'start'])

    # Required VCF columns
    df_vcf = pd.DataFrame({
        'CHROM': df['chr'],
        'POS': df['start'],
        'FILTER': 'PASS',
        'REF': df['reference'],
        'ALT': df['alt']
    })

    # VCF Header
    vcf_header = [
        "#CHROM\tPOS\tFILTER\tREF\tALT"
    ]

    # Group by sample and write VCF per sample
    for sample_id, group in df.groupby("sample"):
        vcf_df = pd.DataFrame({
            'CHROM': group['chr'],
            'POS': group['start'],
            'FILTER': 'PASS',
            'REF': group['reference'],
            'ALT': group['alt']
        })
        
        filename = Path(f"vcf_per_sample/{sample_id}.vcf")
        filename.parent.mkdir(exist_ok=True, parents=True)

        with open(filename, "w") as f:
            for line in vcf_header:
                f.write(line + "\n")
            vcf_df.to_csv(f, sep='\t', index=False, header=False)
    

if __name__ == '__main__':
    file_path = "../data/BRCA_mc3.txt"
    process_dataset(file_path)