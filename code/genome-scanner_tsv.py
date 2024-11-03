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

import importlib
import sys

# Imports a module dynamically and instantiates a class from it
def import_and_instantiate(module_name, class_name, *args, **kwargs):

    module = importlib.import_module(module_name)
    class_ = getattr(module, class_name)
    return class_(*args, **kwargs)


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

    # Import module and instantiate object for specific condition
    c_spondyloarthritis = import_and_instantiate("spondyloarthritis", "Spondyloarthritis", snp_file_name, snp_db_file_name)
    #spondyloarthritis = importlib.import_module("spondyloarthritis")
    #from spondyloarthritis import Spondyloarthritis

    # Create condition spondyloarthritis
    #c_spondyloarthritis = globals()["Spondyloarthritis"](snp_file_name, snp_db_file_name)
    #c_spondyloarthritis = Spondyloarthritis(snp_file_name, snp_db_file_name)

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


