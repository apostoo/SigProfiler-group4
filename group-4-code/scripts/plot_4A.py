from SigProfilerExtractor import sigpro as sig
import SigProfilerExtractor as spe_mod
import shutil 

def run_matrix_96():
    data = sig.importdata("matrix")
    sig.sigProfilerExtractor(
        "matrix",
        "test_matrix_96_output",
        data,
        minimum_signatures=4,
        maximum_signatures=4,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
    )


if __name__ == "__main__":
    run_matrix_96()
    #src_path = "../test_matrix_96_output/SBS96/All_Solutions/SBS96_5_Signatures/Signatures/Signature_plot/SBS_96_plots_S5.pdf"
    src_path = "../test_matrix_96_output/SBS96/All_Solutions/SBS96_4_Signatures/Signatures/Signature_plot/SBS_96_plots_S4.pdf"
    dst_path = "../figures/figure_4A_S4.pdf"
    shutil.copyfile(src_path, dst_path)