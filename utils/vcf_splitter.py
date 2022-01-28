#!/usr/bin/env python
import os
import sys
import argparse
import pandas as pd


def parse_args(argv):
    """ VCF file lift over commands default cli interface """
    help_text = \
        """
        === ISNP VCF Files Formatter ===
        This script takes a main VCF files and splits into sample specific VCF files. 
        """
    args_parser = argparse.ArgumentParser(description=help_text)
    args_parser.add_argument('-i', '--input_vcf', help="<the input vcf file to be split> [mandatory]",
                             action='store', dest="input_vcf", required=True)
    args_parser.add_argument('-o', '--output_dir', help="<output directory for each individual vcf file> [mandatory]",
                             action='store', dest="output_dir", required=True)
    args_parser.add_argument('-p', '--prefix', help="<the prefix to add for chromosome field [optional]>",
                             action='store', dest="prefix", required=False, default="chr")
    return args_parser.parse_args(argv[1:])


def _split_vcf(file, output_dir_path, info_line_id="#CHROM", prefix='chr'):
    """ Split VCF file into individual vcf files """
    with open(file, "r") as read_vcf_file:
        info_lines, headers = [], []
        for line in read_vcf_file:
            if line.startswith(info_line_id):
                headers = line.strip().split('\t')
            elif line.startswith('##'):
                info_lines.append(line.strip())
    assert len(headers) != 0, "Non matching id parameter"
    read_vcf_file_df = pd.read_csv(file, sep='\t', comment='#', header=None)
    read_vcf_file_df.columns = headers
    read_vcf_file_df[info_line_id] = prefix + read_vcf_file_df[info_line_id].astype(str)
    standard_header = headers[0:9]
    sample_ids = headers[9:]
    for sample in sample_ids:
        print(f"Processing: {sample}")
        sample_indices = standard_header + [sample]
        sub_df = read_vcf_file_df[sample_indices]
        sample_file_path = os.path.join(output_dir_path, f"{sample}.vcf")
        with open(sample_file_path, 'w') as out_vcf_file:
            for line in info_lines:
                out_vcf_file.write(line + '\n')
        sub_df.to_csv(sample_file_path, mode='a', sep='\t', index=False, header=True)
    print("Done. :-)")


def split_vcf(argv):
    """ Main logic for vcf splitter """
    args = parse_args(argv)
    _split_vcf(file=args.input_vcf, output_dir_path=args.output_dir, prefix=args.prefix)


if __name__ == "__main__":
    sys.exit(split_vcf(sys.argv))
