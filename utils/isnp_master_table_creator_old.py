import pandas as pd
import argparse
import sys
import os
import re
import logging
import csv

def setup_logging(output_folder):
    """Sets up logging to both console and a file."""
    log_file = "isnp_master_table_creator.log"
    
    # Create logger
    logger = logging.getLogger("ISNP_MasterTable")
    logger.setLevel(logging.INFO)
    
    # Avoid duplicate handlers if setup is called multiple times
    if logger.handlers:
        return logger

    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # File handler
    fh = logging.FileHandler(log_file, mode='w')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger

def parse_args(argv):
    """ Command line interface for the the module """
    help_text = \
        """
        === ISNP Master Table Creator ===
        This script generates a table, which contains all of the TF-target gene and miRNA-target gene interactions per patients.
        """
    parser = argparse.ArgumentParser(description=help_text)

    # Input folder path to learn from
    parser.add_argument("-pf", "--patient_folder",
                        help="<path to the folder where the patient specific folders are> [mandatory]",
                        dest="patient_folder",
                        action="store",
                        required=True)

    results = parser.parse_args(argv)

    return results


def read_in_folder(patient_folder):
    """
    - param folder: patient folder which we use
    - return:

    There should be two different type of interaction.
    # 1: The TF with the uniprot interactions and their target genes
    # 2. The miRNA and their target genes
    """
    logger = setup_logging(patient_folder)
    logger.info(f"Starting master table creation in: {patient_folder}")

    # Data collection
    all_snps_data = []
    all_proteins_data = []
    all_outlines_data = []

    subfolders = [f.path for f in os.scandir(patient_folder) if f.is_dir()]
    logger.info(f"Found {len(subfolders)} patient subfolders.")

    for folder in subfolders:
        patient = os.path.basename(folder)
        filepath = os.path.join(folder, "differences.tsv")
        
        if not os.path.exists(filepath):
            logger.warning(f"File not found: {filepath}. Skipping.")
            continue
            
        logger.info(f"Processing patient: {patient}")
        
        try:
            # Read only necessary columns to save memory and time
            df = pd.read_csv(filepath, sep="\t", header=0)
            
            if df.empty:
                logger.warning(f"Empty differences.tsv for patient {patient}")
                continue

            # Filtering: Keep rows where file != 'mut' OR mutated == True
            # (Matches: if row["file"] == "mut" and row["mutated"] == False: continue)
            df_filtered = df[~((df["file"] == "mut") & (df["mutated"] == False))].copy()

            if df_filtered.empty:
                continue

            # SNPs
            snps = df_filtered["snp"].astype(str).str.strip().unique()
            for snp in snps:
                all_snps_data.append((snp, patient))

            # Proteins (target)
            proteins = df_filtered["target"].astype(str).str.strip().unique()
            for protein in proteins:
                all_proteins_data.append((protein, patient))

            # Outlines
            # Vectorized outline construction
            df_filtered["outline_str"] = (
                df_filtered["source"].astype(str).str.strip() + "\t" +
                df_filtered["target"].astype(str).str.strip() + "\t" +
                df_filtered["snp"].astype(str).str.strip() + "\t" +
                df_filtered["file"].astype(str).str.strip() + "\t" +
                df_filtered["tool"].astype(str).str.strip()
            )
            outlines = df_filtered["outline_str"].unique()
            for outline in outlines:
                all_outlines_data.append((outline, patient))

        except Exception as e:
            logger.error(f"Error processing {patient}: {e}")

    logger.info("Building final dataframes...")

    def create_and_format_df(data, index_name):
        if not data:
            return pd.DataFrame(columns=["ID"]).set_index("ID")
        
        df_long = pd.DataFrame(data, columns=[index_name, "Patient"])
        df_long["Value"] = 1
        
        # Pivot to get Patient columns and Item index
        df_pivot = df_long.pivot_table(index=index_name, columns="Patient", values="Value", aggfunc="max").fillna(0).astype(int)
        
        # Upper case index as in original script
        df_pivot.index = df_pivot.index.str.upper()
        
        return df_pivot

    # Create DataFrames
    rs_df = create_and_format_df(all_snps_data, "SNP")
    proteins_df = create_and_format_df(all_proteins_data, "ID")
    outline_df = create_and_format_df(all_outlines_data, "Source\tTarget\tSNP\tMutated\tInteraction_source")

    # Rename axis as in original (though pivot already sets the index name)
    rs_df = rs_df.rename_axis("SNP")
    proteins_df = proteins_df.rename_axis("ID")
    outline_df = outline_df.rename_axis("Source\tTarget\tSNP\tMutated\tInteraction_source")

    logger.info(f"Final matrix sizes: SNPs: {rs_df.shape}, Proteins: {proteins_df.shape}, Outlines: {outline_df.shape}")

    # Write out the results
    logger.info("Writing results to TSV files...")
    
    # Original used quoting=False which is 0 (minimal). Using csv.QUOTE_MINIMAL for clarity.
    proteins_df.to_csv(os.path.join(patient_folder, "affected_proteins.tsv"), sep="\t", quoting=csv.QUOTE_MINIMAL)
    outline_df.to_csv(os.path.join(patient_folder, "affected_proteins_TFs_mirs.tsv"), sep="\t", quoting=csv.QUOTE_MINIMAL)
    rs_df.to_csv(os.path.join(patient_folder, "SNPs.tsv"), sep="\t", quoting=csv.QUOTE_MINIMAL)

    logger.info("Master table creation completed successfully.")

def main(argv):
    """Main function for the script"""
    args = parse_args(argv)
    read_in_folder(args.patient_folder)

if __name__ == "__main__":
    main(sys.argv[1:])

