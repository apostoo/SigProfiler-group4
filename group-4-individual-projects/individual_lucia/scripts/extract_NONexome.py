from SigProfilerExtractor import sigpro as sig
import csv
import numpy as np
import matplotlib.pyplot as plt

def extract_wng(data,n):
        sig.sigProfilerExtractor(
        "vcf",
        "individual_lucia/wNg"+str(n),
        data,
        minimum_signatures=n,
        maximum_signatures=n,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
        )


"""
    Run the non-exome for multiple value of n to find best number of signatures on non-exome data
"""
def run_for_all_signature_numbers(data):

    # Iterate up until 15 signatures
    for n in range(1, 16):

        extract_wng(data, n)


def build_csv():
    output_file = 'individual_lucia/aggregated_NONexome_solutions.csv'
    with open(output_file, 'a',newline='') as final_csv:
        writer = csv.writer(final_csv)

        # Save header
        with open("individual_lucia/wNg1/SBS96/All_solutions_stat.csv", 'r') as output_csv:
            reader = csv.reader(output_csv)
            writer.writerow(next(reader))

        # Look up solution statistics
        for i in range(1, 16):
            path = "individual_lucia/wNg"+str(i)+"/SBS96/All_solutions_stat.csv"
            # print(path)
            with open(path, 'r') as output_csv:
                reader = list(csv.reader(output_csv))
                if len(reader) > 1:
                    writer.writerow(reader[-1])

    return output_file      

"""
    Get the signatures number, stability and frobenius% from the csv.
"""
def get_csv_values(file):
    signatures = []
    stability = []
    frobenius = []

    with open(file, 'r') as data:
        reader = csv.reader(data)
        next(reader)

        for row in reader:

            # Might contain *
            raw_sig = row[0].strip()
            if raw_sig.endswith('*'):
                raw_sig = raw_sig[:-1]

            signatures.append(int(raw_sig))
            stability.append(float(row[1]))
            frobenius.append(float(row[5][:-1])/100)
        
    return (signatures, stability, frobenius)

"""
    Plot the lines according to the 3 variables, just like in the SigProfiler paper.
"""
def plot_lines(signatures, reproducibility, frobenius):
    fig, ax1 = plt.subplots()

    # Plot x-axis for signatures number
    ax1.set_xlabel('Number of mutational signatures')
    ax1.set_xticks(np.arange(1, 16, 1.0))

    # Plot line corresponding to signature reproducibility.
    color = 'tab:red'
    ax1.set_ylabel('Signature reproducibility', color=color)
    ax1.plot(signatures, reproducibility, color=color, marker='o', label='Signature reproducibility')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_ylim(0.0, 1.05)
    ax1.set_xlim(1, 15)

    # Plot line corresponding to Frobenius reconstruction error.
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Frobenius reconstruction error', color=color)
    ax2.plot(signatures, frobenius, color=color, marker='s', label='Frobenius reconstruction error')
    ax2.tick_params(axis='y', labelcolor=color)

    # plt.title('Mutational Signature Analysis')
    plt.tight_layout()

    plt.savefig("individual_lucia/NONExome_convergence.pdf")


if __name__ == "__main__":

    vcf_data =  "individual_lucia/non_exome_data/"

    # Find N
    # signatures_csv = run_for_all_signature_numbers(vcf_data)
    signatures_csv = 'individual_lucia/aggregated_NONexome_solutions.csv'

    # Combine solutions
    # build_csv()

    # Get necessary values from the csv 
    signatures, reproducibility, frobenius = get_csv_values(signatures_csv)

    # Plot the results
    plot_lines(signatures, reproducibility, frobenius)