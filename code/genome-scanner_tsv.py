#!/usr/bin/env python3
"""
Scans TSV files for SNPs from a database and interprets results
"""
__author__ = "Melanie Senn"
__copyright__ = "Copyright 2024"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Melanie Senn"
__email__ = "melanie.senn@gmail.com"

from spondyloarthritis import Spondyloarthritis
import sys

if __name__ == '__main__':

     # Get command line arguments
    args = sys.argv

    # File names
    snp_file_name = None
    snp_db_file_name = "../db/snpdb_spondyloarthritis.csv"

    # Check if there are command line arguments
    if len(args) > 1:
        # Get the first argument as the user's SNP filename
        snp_file_name = args[1]

    # Create condition spondyloarthritis
    c_spondyloarthritis = Spondyloarthritis(snp_file_name, snp_db_file_name)

    # Get content from TSV files
    c_spondyloarthritis.get_tsv_content()

    # Load SNP database
    c_spondyloarthritis.get_snp_db()

    # Get results for SNPs found in TSV content and SNP database
    c_spondyloarthritis.get_snp_results()

    # Summarize results
    c_spondyloarthritis.summarize_results()

    # Save results (csv and json)
    c_spondyloarthritis.save_results("results.csv", "results.json")

    print('\nScan finished')


