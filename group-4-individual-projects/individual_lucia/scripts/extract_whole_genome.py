from SigProfilerExtractor import sigpro as sig

# Extracts 4 signatures as indicated by previous group work
def extract_wsg(data):
        sig.sigProfilerExtractor(
        "vcf",
        "individual_lucia/wSg_4",
        data,
        minimum_signatures=1,
        maximum_signatures=4,
        nmf_replicates=5,
        min_nmf_iterations=100,
        max_nmf_iterations=1000,
        nmf_test_conv=100,
        )

        
        
if __name__ == "__main__":
        vcf_data =  "datasets/unzipped/21BRCA/21BRCA/21BRCA_vcf/"
        extract_wsg(vcf_data)