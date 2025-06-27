import pandas as pd


     
if __name__ == "__main__":

    # WES 
    df_e = pd.read_csv('individual_lucia\wEg2\SBS96\Suggested_Solution\SBS96_De-Novo_Solution\Signatures\SBS96_De-Novo_Signatures.txt', sep='\t')
   
    # WSG data
    df_s = pd.read_csv('individual_lucia\wSg_4\SBS96\Suggested_Solution\SBS96_De-Novo_Solution\Signatures\SBS96_De-Novo_Signatures.txt', sep='\t')


    # WNS data
    df = pd.read_csv('individual_lucia\wNg4\SBS96\Suggested_Solution\SBS96_De-Novo_Solution\Signatures\SBS96_De-Novo_Signatures.txt', sep='\t')
   
   # Example difference
    diff = (df_e['SBS96A'] -df_s['SBS96A']).values

    print(diff)


    