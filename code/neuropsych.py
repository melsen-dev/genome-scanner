#!/usr/bin/env python3
"""
Describes condition Neuropsychatric Disorders
"""
__author__ = "Melanie Senn"
__copyright__ = "Copyright 2025"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Melanie Senn"
__email__ = "melanie.senn@gmail.com"

from condition import *
import json

# Name constants
ASSOCIATIONS_DIAGNOSIS = ['Average risk', 'Small increase in risk', 'Increased risk']
V_ASD = 'Autism Spectrum Disorder'
V_TS = 'Tourettes Syndrome'
V_BIP = 'Bipolar Disorder'
V_SCZ = 'Schizophrenia'
V_MDD = 'Major Depressive Disorder'
V_ANO = 'Anorexia Nervosa'

class Neuropsych(Condition):
    def __init__(self, snp_file_name, snp_db_file_name):
        super().__init__(snp_file_name, snp_db_file_name)

        self.diagnostic_risk_association = []

    # Get risk association for result row
    def get_risk_association(self, result_row):
        association = ''
        risk = 0
        application = result_row[C_APPLICATION]
    
        # Risk = 0: if risk allele not found in genotype
        # Risk = 1: if risk allele found once in genotype
        # Risk = 2: if risk allele found twice in genotype
        risk = result_row[C_GENOTYPE].count(result_row[C_RISK_ALLELE])
        max_risk = 2
        assert((risk < 2) or (risk > 0))
        association = ASSOCIATIONS_DIAGNOSIS[risk]
        
        return risk, max_risk, association
    
    # Summarize results from single SNPs into conclusion
    def summarize_results(self):

        # For each disorder
        disorders = [V_ASD, V_TS, V_BIP, V_SCZ, V_MDD, V_ANO]
        for disorder in disorders:
            self.diagnostic_risk_association.append("Unknown risk")
        for disorder_count, disorder in enumerate(disorders):

            # Diagnostic summary
            diagnostic_results = self.snp_results.loc[self.snp_results[C_APPLICATION] == disorder]
            if len(diagnostic_results):
                diagnostic_risk = diagnostic_results[C_RISK].sum()
                diagnostic_max_risk = diagnostic_results[C_MAX_RISK].sum()
                diagnostic_rel_risk = diagnostic_risk / diagnostic_max_risk

                # Average risk
                if diagnostic_rel_risk == 0.0:
                    diagnostic_risk_association = ASSOCIATIONS_DIAGNOSIS[0]

                # Small increase in risk
                elif diagnostic_rel_risk <= 0.5:
                    diagnostic_risk_association = ASSOCIATIONS_DIAGNOSIS[1]

                # Increased risk
                else:
                    diagnostic_risk_association = ASSOCIATIONS_DIAGNOSIS[2]

                self.diagnostic_risk_association[disorder_count] = diagnostic_risk_association

            # Summarize results from single SNPs for conclusion
            print('\nYour diagnostic insights for ' + disorder + ': ' + self.diagnostic_risk_association[disorder_count])

    # Save results (csv and json)
    def save_results(self, results_csv_file_name, results_json_file_name):

        # csv file
        results = self.snp_results.drop(columns=[C_CONDITION, C_RISK, C_MAX_RISK])
        print('\nResults\n' + results.to_markdown())
        results.to_csv(results_csv_file_name, index=False, sep='\t')

        # json file
        result_dic = {
            "dianogstic_score_asd": self.diagnostic_risk_association[0],
            "dianogstic_score_ts": self.diagnostic_risk_association[1],
            "dianogstic_score_bip": self.diagnostic_risk_association[2],
            "dianogstic_score_scz": self.diagnostic_risk_association[3],
            "dianogstic_score_mdd": self.diagnostic_risk_association[4],
            "dianogstic_score_ano": self.diagnostic_risk_association[5],
            "table_csv" : self.snp_results.to_csv(None, index=False, sep='\t')
        }

        # Convert and write JSON object to file
        with open(results_json_file_name, "w") as outfile:
            json.dump(result_dic, outfile)
