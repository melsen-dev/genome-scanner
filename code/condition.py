#!/usr/bin/env python3
"""
Describes a condition (medical, wellness etc.)
"""
__author__ = "Melanie Senn"
__copyright__ = "Copyright 2024"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Melanie Senn"
__email__ = "melanie.senn@gmail.com"

from abc import abstractmethod
import pandas as pd
from zipfile import ZipFile
from tkinter import filedialog as fd

# Name constants
C_APPLICATION = 'Application'
C_CONDITION = 'Condition'
C_SNP = 'SNP'
C_GENE = 'Gene/locus'
C_GENOTYPE = 'Genotype'
C_RISK_ALLELE = 'Risk allele'
C_PROTECTIVE_ALLELE = 'Protective allele'
C_REFERENCE = 'Reference'
C_ASSOCIATION = 'Association'
C_RISK = 'Risk'
C_MAX_RISK = 'Max risk'
V_DIAGNOSIS = 'Diagnosis'
V_TREATMENT = 'Treatment'

class Condition:
    def __init__(self, snp_file_name, snp_db_file_name):
        self.snp_file_name = snp_file_name
        self.snp_db_file_name = snp_db_file_name
        self.tsv_content = None
        self.snp_db = None
        self.snp_results = None

    # Get TSV files content from user selection
    def get_tsv_content(self):

        if self.snp_file_name is None:
            file_types = (('zip files', '*.zip'), ('All files', '*.*'))
            self.snp_file_name = fd.askopenfilename(title='Open a ZIP file', initialdir='/', filetypes=file_types)

        # Check if file is a zip file
        if self.snp_file_name.endswith('.zip'):
            
            # Get TSV content from zip file
            zip_file = ZipFile(self.snp_file_name)
            tsv_content = pd.concat([pd.read_csv(zip_file.open(i.filename), sep='\t', header=None) for i in zip_file.infolist() if i.compress_size > 0], ignore_index=True)
        else:
            tsv_content = pd.read_csv(self.snp_file_name, sep='\t', header=None)

        self.tsv_content = tsv_content
    
    # Load SNP database
    def get_snp_db(self):
        self.snp_db = pd.read_csv(self.snp_db_file_name)

    # Get results from TSV files for database SNPs
    def get_snp_results(self):

        # Create empty results table
        result_columns=[C_CONDITION, C_APPLICATION, C_SNP, C_GENE, C_GENOTYPE, C_RISK_ALLELE, C_PROTECTIVE_ALLELE, C_ASSOCIATION, C_REFERENCE]
        results = pd.DataFrame(columns=result_columns)

        # Search tsv data for SNPs
        for db_index in self.snp_db.index:
            current_SNP = self.snp_db[C_SNP][db_index]
            print(f"Scanning for SNP {current_SNP} ({db_index + 1} out of {len(self.snp_db.index)})")
            SNP_match = self.tsv_content[0] == current_SNP
            if SNP_match.sum() == 1: # check for unique match
                print(f"Found SNP {current_SNP}")
                row = self.tsv_content.loc[SNP_match]
                result_entry = len(results.index)
                results.loc[result_entry, C_SNP] = current_SNP
                results.loc[result_entry, C_CONDITION] = self.snp_db[C_CONDITION][db_index]
                results.loc[result_entry, C_APPLICATION] = self.snp_db[C_APPLICATION][db_index]
                results.loc[result_entry, C_GENE] = self.snp_db[C_GENE][db_index]
                results.loc[result_entry, C_REFERENCE] = self.snp_db[C_REFERENCE][db_index]
                results.loc[result_entry, C_RISK_ALLELE] = self.snp_db[C_RISK_ALLELE][db_index]
                results.loc[result_entry, C_PROTECTIVE_ALLELE] = self.snp_db[C_PROTECTIVE_ALLELE][db_index]
                results.loc[result_entry, C_GENOTYPE] = row[3].values[0] # genotype
                risk, max_risk, association = self.get_risk_association(results.loc[result_entry]) # risk association
                results.loc[result_entry, C_RISK] = risk
                results.loc[result_entry, C_MAX_RISK] = max_risk
                results.loc[result_entry, C_ASSOCIATION] = association

        self.snp_results = results

    # Get risk association for result row
    @abstractmethod
    def get_risk_association(self, result_row):
        pass

    # Summarize results from single SNPs into conclusion
    @abstractmethod
    def summarize_results(self):
        pass

    # Save results (csv and json)
    def save_results(self, results_csv_file_name, results_json_file_name):

        # csv file
        results = self.snp_results.drop(columns=[C_APPLICATION, C_CONDITION, C_RISK, C_MAX_RISK])
        print('\nResults\n' + results.to_markdown())
        results.to_csv(results_csv_file_name, index=False, sep='\t')
