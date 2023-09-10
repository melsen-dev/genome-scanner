#!/usr/bin/env python3
"""
Scans VCF file for SNPs from a database and interprets results
"""
__author__ = "Melanie Senn"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Melanie Senn"
__email__ = "melanie.senn@gmail.com"

from pysam import VariantFile
import numpy as np
import pandas as pd
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

# Get genotype for record
def get_genotype(record):
    
    # Currently do not handle delete/insert
    # Reference: reference genome
    # Alt: studied sample
    # GT: 0 reference allele, 1 alternative allele
    # If ref_len > alt_len --> delete, if ref_len < alt_len --> insert
    ref_len = len(record.ref)   
    alt_len = len(record.alts)
    assert(ref_len == alt_len)  

    # For details on genotype see pysam documentation
    sample_id = record.samples.keys()[0]
    genotype_binary = record.samples[sample_id]['GT']
    genotype = 'II' # Invalid genotype
    assert(genotype_binary != (0, 0)) # (0, 0) should never happen
    if (genotype_binary == (1, 1)): # Homozygous alt
        genotype = record.alts[0] + record.alts[0]
    elif (genotype_binary == (0, 1)): # Heterozygous
        genotype = record.ref + record.alts[0]
    elif (genotype_binary == (0, 0)): # Homozygous ref
        genotype = record.ref + record.ref
    #print(f"Binary genotype {genotype_binary}")
    
    return genotype

# Get results from VCF file for database SNPs
def get_snp_results(database, vcf_file):
    
    # Create empty results table
    result_columns=[C_CONDITION, C_APPLICATION, C_SNP, C_GENE, C_GENOTYPE, C_RISK_ALLELE, C_PROTECTIVE_ALLELE, C_ASSOCIATION, C_REFERENCE]
    results = pd.DataFrame(columns=result_columns)

    # Find results in VCF file for each SNP in database 
    for db_index in database.index:
        current_SNP = database[C_SNP][db_index]
        print(f"Scanning for SNP {current_SNP} ({db_index + 1} out of {len(database.index)})")
        records = vcf_file.fetch()
        for record in records:
            if (record.id == current_SNP):
                print(f"Found SNP {current_SNP}")
                result_entry = len(results.index)
                results.loc[result_entry, C_SNP] = current_SNP
                results.loc[result_entry, C_CONDITION] = database[C_CONDITION][db_index]
                results.loc[result_entry, C_APPLICATION] = database[C_APPLICATION][db_index]
                results.loc[result_entry, C_GENE] = database[C_GENE][db_index]
                results.loc[result_entry, C_REFERENCE] = database[C_REFERENCE][db_index]
                results.loc[result_entry, C_RISK_ALLELE] = database[C_RISK_ALLELE][db_index]
                results.loc[result_entry, C_PROTECTIVE_ALLELE] = database[C_PROTECTIVE_ALLELE][db_index]
                results.loc[result_entry, C_GENOTYPE] = get_genotype(record)
                results.loc[result_entry, C_ASSOCIATION] = get_association(results.loc[result_entry])
                break

    return results

# Load SNP database
def get_snpdb(db_filename):

    database = pd.read_csv(db_filename)
    return database

# Get VCF file from user selection
def get_vcf_file():

    file_types = (('vcf files', '*.vcf'), ('All files', '*.*'))
    file_name = fd.askopenfilename(title='Open a VCF file', initialdir='/', filetypes=file_types)
    file = VariantFile(file_name)

    return file

if __name__ == '__main__':
    
    # Get VCF file
    vcf_file = get_vcf_file()

    # Load SNP database
    snpdb = get_snpdb('../db/snpdb_sa.csv')

    # Get results for SNPs found in VCF file and SNP database
    results = get_snp_results(snpdb, vcf_file)

    # Print results
    print('\nResults\n' + results.to_markdown())

    # Print SNP database
    print('\nSNP Database\n' + snpdb.to_markdown())

    print('\nScan finished')


