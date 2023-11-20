#!/usr/bin/env python3
"""
Scans TSV files for SNPs from a database and interprets results
"""
__author__ = "Melanie Senn"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Melanie Senn"
__email__ = "melanie.senn@gmail.com"

import numpy as np
import pandas as pd
from tkinter import filedialog as fd
import glob
from zipfile import ZipFile

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
V_DIAGNOSIS = 'Diagnosis'
V_TREATMENT = 'Treatment'
ASSOCIATIONS_DIAGNOSIS = ['Average risk', 'Small increase in risk', 'Increased risk']
ASSOCIATIONS_TREATMENT = ['clinical response']

# Get association for result row
def get_association(result_row):

    association = ''
    application = result_row[C_APPLICATION]
    if application == V_DIAGNOSIS:
        
        # Risk = 0: if risk allele not found in genotype
        # Risk = 1: if risk allele found once in genotype
        # Risk = 2: if risk allele found twice in genotype
        risk = result_row[C_GENOTYPE].count(result_row[C_RISK_ALLELE])
        assert((risk < 2) or (risk > 0))
        association = ASSOCIATIONS_DIAGNOSIS[risk]
    elif V_TREATMENT in application:
        
        # Break down application string
        application, medication, response = application.split('_')
        association = 'Unknown ' + ASSOCIATIONS_TREATMENT[0] + ' for ' + medication
        risk_allele = result_row[C_RISK_ALLELE]
        len_risk_allele = len(risk_allele)
        genotype = result_row[C_GENOTYPE]
        
        # One risk allele given by snpdb (paper reference)
        if len_risk_allele == 1:     
            risk = genotype.count(risk_allele)
            assert((risk < 2) or (risk > 0))
            if risk > 0:
                association = response + ' ' + ASSOCIATIONS_TREATMENT[0] + ' for ' + medication
        
        # Two combined risk alleles as genotype given by snpdb (paper reference)
        elif len_risk_allele == 2:  
            reversed_genotype = genotype[::-1]
            if (risk_allele == genotype) or (risk_allele == reversed_genotype):
                association = response + ' ' + ASSOCIATIONS_TREATMENT[0] + ' for ' + medication
    
    return association

# Get results from TSV files for database SNPs
def get_snp_results(database, tsv_files):
    
    # Create empty results table
    result_columns=[C_CONDITION, C_APPLICATION, C_SNP, C_GENE, C_GENOTYPE, C_RISK_ALLELE, C_PROTECTIVE_ALLELE, C_ASSOCIATION, C_REFERENCE]
    results = pd.DataFrame(columns=result_columns)

    # Search tsv data for SNPs
    for db_index in database.index:
        current_SNP = database[C_SNP][db_index]
        print(f"Scanning for SNP {current_SNP} ({db_index + 1} out of {len(database.index)})")
        SNP_match = tsv_content[0] == current_SNP
        if SNP_match.sum() == 1: # check for unique match
            print(f"Found SNP {current_SNP}")
            row = tsv_content.loc[SNP_match]
            result_entry = len(results.index)
            results.loc[result_entry, C_SNP] = current_SNP
            results.loc[result_entry, C_CONDITION] = database[C_CONDITION][db_index]
            results.loc[result_entry, C_APPLICATION] = database[C_APPLICATION][db_index]
            results.loc[result_entry, C_GENE] = database[C_GENE][db_index]
            results.loc[result_entry, C_REFERENCE] = database[C_REFERENCE][db_index]
            results.loc[result_entry, C_RISK_ALLELE] = database[C_RISK_ALLELE][db_index]
            results.loc[result_entry, C_PROTECTIVE_ALLELE] = database[C_PROTECTIVE_ALLELE][db_index]
            results.loc[result_entry, C_GENOTYPE] = row[3].values[0] # genotype
            results.loc[result_entry, C_ASSOCIATION] = get_association(results.loc[result_entry]) # association

    return results

# Load SNP database
def get_snpdb(db_filename):

    database = pd.read_csv(db_filename)
    return database

# Get TSV files content from user selection
def get_tsv_content():

    file_types = (('zip files', '*.zip'), ('All files', '*.*'))
    file_name = fd.askopenfilename(title='Open a ZIP file', initialdir='/', filetypes=file_types)

    # Get TSV content from zip file
    zip_file = ZipFile(file_name)
    tsv_content = pd.concat([pd.read_csv(zip_file.open(i.filename), sep='\t', header=None) for i in zip_file.infolist() if i.compress_size > 0], ignore_index=True)
    return tsv_content

if __name__ == '__main__':
    
    # Get content from TSV files
    tsv_content = get_tsv_content()

    # Load SNP database
    snpdb = get_snpdb('../db/snpdb_sa.csv')

    # Get results for SNPs found in TSV content and SNP database
    results = get_snp_results(snpdb, tsv_content)

    # Print SNP database
    #print('\nSNP Database\n' + snpdb.to_markdown())

    # Print results
    results = results.drop(columns=[C_APPLICATION, C_CONDITION])
    print('\nResults\n' + results.to_markdown())
    results.to_csv('results.csv', index=False, sep='\t')

    print('\nScan finished')


