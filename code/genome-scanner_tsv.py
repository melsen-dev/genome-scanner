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
import json

# Import a module dynamically and instantiates an object
def import_and_instantiate(module_name, class_name, *args, **kwargs):

    module = importlib.import_module(module_name)
    class_ = getattr(module, class_name)
    return class_(*args, **kwargs)

if __name__ == '__main__':

    # SNP file name
    snp_file_name = None

    # Check command line arguments
    args = sys.argv
    if len(args) > 1:
        # Get the first argument as the condition's config
        condition_config = args[1]
    else:
        print("You need to provide at least the config file as the first command line parameter")
        exit(0)
    if len(args) > 2:
        # Get the second argument as the user's SNP file
        snp_file_name = args[2]

    with open(condition_config, "r") as config_file:
        config_data = json.load(config_file)

    # Import module and instantiate object for specific condition with config data
    condition = import_and_instantiate(config_data["python_module"], config_data["class_constructor"],
                                       snp_file_name, config_data["snp_database"])

    # Get content from TSV files
    condition.get_tsv_content()

    # Load SNP database
    condition.get_snp_db()

    # Get results for SNPs found in TSV content and SNP database
    condition.get_snp_results()

    # Summarize results
    condition.summarize_results()

    # Save results (csv and json)
    condition.save_results("results.csv", "results.json")

    print("\nScan finished")
