import pandas as pd
import numpy as np
import os
import glob

def calculate_tpm(count_file):
    """Calculates TPM from samtools idxstats output."""
    # idxstats: col0=ID, col1=Length, col2=MappedReads
    df = pd.read_csv(count_file, sep='\t', header=None, 
                     names=['GeneID', 'Length', 'Counts', 'Unmapped'])
    df = df[df['GeneID'] != '*']  # Remove unmapped row
    
    # RPK = Counts / (Length in kb)
    rpk = df['Counts'] / (df['Length'] / 1000.0)
    # TPM = RPK / (Sum of RPK / 10^6)
    scaling_factor = rpk.sum() / 1000000.0
    if scaling_factor == 0:
        df['TPM'] = 0
    else:
        df['TPM'] = rpk / scaling_factor
    return df[['GeneID', 'TPM']]

def main():
    print("Step 1: Processing abundance files...")
    count_files = glob.glob('3-gene_catalog/*.counts')
    abundance_list = []
    
    for f in count_files:
        sample_name = os.path.basename(f).replace('.counts', '')
        tpm_df = calculate_tpm(f)
        tpm_df.columns = ['GeneID', sample_name]
        abundance_list.append(tpm_df.set_index('GeneID'))
    
    # Create the TPM Matrix
    abundance_matrix = pd.concat(abundance_list, axis=1).fillna(0)

    print("Step 2: Loading annotation files...")
    # Load KEGG (Assuming diamond blastp -f 6 format)
    # Adjust names/paths based on your actual filenames
    kegg_df = pd.read_csv('4-tax_kegg/kegg.txt', sep='\t', header=None,
                          usecols=[0, 1], names=['GeneID', 'KEGG_KO']).drop_duplicates('GeneID')

    # Load Taxonomy (Assuming NR/CAT results)
    tax_df = pd.read_csv('4-tax_kegg/nr.daa_tax.txt', sep='\t', 
                         header=None, names=['GeneID', 'Taxonomy']).drop_duplicates('GeneID')

    print("Step 3: Merging all data into a Master Matrix...")
    # Merge Annotation first
    annotations = pd.merge(kegg_df, tax_df, on='GeneID', how='outer')
    
    # Merge with Abundance Matrix
    master_matrix = pd.merge(annotations, abundance_matrix, left_on='GeneID', right_index=True, how='right')

    # Save to CSV
    output_file = 'final_gene_abundance_metadata_matrix.csv'
    master_matrix.to_csv(output_file)
    print(f"Success! Master matrix saved to {output_file}")

if __name__ == "__main__":
    main()