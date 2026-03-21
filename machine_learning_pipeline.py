import pandas as pd
import numpy as np
# ---------- STEP 0: import ClinVar (NCBI dataset) and exploring the dataset: -------
#var_summarized=pd.read_csv("~/Desktop/colon_cancer_metaG/variant_summary.txt", sep = "\t", dtype={"Chromosome": str})

# 0. Subset the first 2k samples (1k head, 1k tail, then combine)
#var_df=pd.concat([var_summarized.head(1000),
#                 var_summarized.tail(1000)], ignore_index = True)
#var_df.to_csv("~/Desktop/colon_cancer_metaG/variant_subset.csv", index = False)

var_summarized=pd.read_csv("~/Desktop/colon_cancer_metaG/variant_subset.csv")
# 1. What does one row look like?
var_summarized.iloc[0]

# 2. What are the unique clinical significance labels
var_summarized["ClinicalSignificance"].unique()

# 3. How many rows have "Benign" or "Pathogenic"?
print(var_summarized[var_summarized["ClinicalSignificance"] == "Benign"].value_counts())
print(var_summarized[var_summarized["ClinicalSignificance"] == "Pathogenic"].value_counts())

# 4. How big is the dataset? - number of rows vs columns
var_summarized.shape #(8907105 rows, 43 columns)

# 5. What information is available? - 5 columns that are important to the project
#important columns:GeneID, Type, ClinicalSignificance, PhenotypeList, OriginSimple, Chromosome
var_summarized_subset=var_summarized[["GeneID", "SCVsForAggregateOncogenicityClassification", "ClinicalSignificance", "PhenotypeList", "OriginSimple", "PositionVCF"]]

# 6. What are we predicting?- predictor vars from the FULL dataset:
outcome_var=var_summarized["ClinicalSignificance"]

# 7. How clean is the data?
#print(var_summarized.isna().values.any()) #check for NA - no NA values
var_summarized.isnull().any #no NA values

# 8. What types of variants exist?
var_summarized["Type"].unique()

# 9. How imbalanced is your target?
var_summarized["Type"].value_counts() #snv is very dominant and complex is only 1

# 10.  Where are the variants located?
print(var_summarized.columns)
print(var_summarized["OriginSimple"].unique()) #'germline', 'not applicable', 'germline/somatic', 'unknown'

