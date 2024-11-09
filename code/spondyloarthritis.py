#!/usr/bin/env python3
"""
Describes condition Spondyloarthritis
"""
__author__ = "Melanie Senn"
__copyright__ = "Copyright 2024"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Melanie Senn"
__email__ = "melanie.senn@gmail.com"

from condition import *
import json

# Name constants
V_TREATMENT_TNF_POS = 'Treatment_TNF-Inhibitor_Positive'
V_TREATMENT_MTX_POS = 'Treatment_Methotrexate_Positive'
V_TREATMENT_CLIN_RESP_POS = ' for positive clinical response'
V_TREATMENT_CLIN_RESP = 'clinical response'
ASSOCIATIONS_DIAGNOSIS = ['Average risk', 'Small increase in risk', 'Increased risk']
ASSOCIATIONS_TREATMENT_OPPORTUNITY = ['Low opportunity', 'High opportunity']

class Spondyloarthritis(Condition):
    def __init__(self, snp_file_name, snp_db_file_name):
        super().__init__(snp_file_name, snp_db_file_name)
        
        self.diagnostic_risk_association = None
        self.treatment_risk_association_tnf_pos = None
        self.treatment_risk_association_mtx_pos = None

    # Get risk association for result row
    def get_risk_association(self, result_row):
        association = ''
        risk = 0
        application = result_row[C_APPLICATION]
        
        if application == V_DIAGNOSIS:
            # Risk = 0: if risk allele not found in genotype
            # Risk = 1: if risk allele found once in genotype
            # Risk = 2: if risk allele found twice in genotype
            risk = result_row[C_GENOTYPE].count(result_row[C_RISK_ALLELE])
            max_risk = 2
            assert((risk < 2) or (risk > 0))
            association = ASSOCIATIONS_DIAGNOSIS[risk]
        
        elif V_TREATMENT in application:

            # Break down application string
            application, medication, response = application.split('_')
            association = 'Unknown ' + V_TREATMENT_CLIN_RESP + ' for ' + medication
            risk_allele = result_row[C_RISK_ALLELE]
            len_risk_allele = len(risk_allele)
            genotype = result_row[C_GENOTYPE]

            # One risk allele given by snpdb (paper reference)
            if len_risk_allele == 1:
                risk = genotype.count(risk_allele)
                max_risk = 2
                assert((risk < 2) or (risk > 0))
                if risk > 0:
                    association = response + ' ' + V_TREATMENT_CLIN_RESP + ' for ' + medication

            # Two combined risk alleles as genotype given by snpdb (paper reference)
            elif len_risk_allele == 2:
                reversed_genotype = genotype[::-1]
                max_risk = 1
                if (risk_allele == genotype) or (risk_allele == reversed_genotype):
                    association = response + ' ' + V_TREATMENT_CLIN_RESP + ' for ' + medication
                    risk = 1

        return risk, max_risk, association
    
    # Summarize results from single SNPs into conclusion
    def summarize_results(self):

        # Diagnostic summary
        diagnostic_results = self.snp_results.loc[self.snp_results[C_APPLICATION] == V_DIAGNOSIS]
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

        # Treatment summary
        # Get relative risk for treatment options
        def get_rel_risk(treatment_application):
            treatment_results = self.snp_results.loc[self.snp_results[C_APPLICATION] == treatment_application]
            treatment_risk = treatment_results[C_RISK].sum()
            treatment_max_risk = treatment_results[C_MAX_RISK].sum()
            treatment_rel_risk = treatment_risk / treatment_max_risk

            # Low treatment opportunity
            if treatment_rel_risk <= 0.3:
                treatment_risk_association = ASSOCIATIONS_TREATMENT_OPPORTUNITY[0]

            # High treatment opportunity
            else:
                treatment_risk_association = ASSOCIATIONS_TREATMENT_OPPORTUNITY[1]

            return treatment_risk_association

        # TNF Inhibitor
        treatment_risk_association_tnf_pos = get_rel_risk(V_TREATMENT_TNF_POS)

        # Methotrexate
        treatment_risk_association_mtx_pos = get_rel_risk(V_TREATMENT_MTX_POS)

        self.diagnostic_risk_association = diagnostic_risk_association
        self.treatment_risk_association_tnf_pos = treatment_risk_association_tnf_pos
        self.treatment_risk_association_mtx_pos = treatment_risk_association_mtx_pos
    
        # Summarize results from single SNPs for conclusion
        print('\nYour diagnostic insights: ' + diagnostic_risk_association)
        print('\nYour TNF inhibitor treatment insights: ' + treatment_risk_association_tnf_pos + V_TREATMENT_CLIN_RESP_POS)
        print('\nYour Methotrexate treatment insights: ' + treatment_risk_association_mtx_pos + V_TREATMENT_CLIN_RESP_POS)

    # Save results (csv and json)
    def save_results(self, results_csv_file_name, results_json_file_name):

        super().save_results(results_csv_file_name, results_json_file_name)

        # json file
        result_dic = {
            "dianogstic_score": self.diagnostic_risk_association,
            "tnf_treatment_score": self.treatment_risk_association_tnf_pos + V_TREATMENT_CLIN_RESP_POS,
            "mr_treatment_score" : self.treatment_risk_association_mtx_pos + V_TREATMENT_CLIN_RESP_POS,
            "table_csv" : self.snp_results.to_csv(None, index=False, sep='\t')
        }

        # Convert and write JSON object to file
        with open(results_json_file_name, "w") as outfile:
            json.dump(result_dic, outfile)
