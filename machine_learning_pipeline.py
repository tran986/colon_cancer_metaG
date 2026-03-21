import pandas as pd
# ---------- STEP 1: import ClinVar (NCBI dataset) and exploring the dataset: -------
var_summarized=pd.read_csv("~/Desktop/colon_cancer_metaG/variant_summary.txt.gz", sep = "\t", dtype={"Chromosome": str})
print(var_summarized.head(5))

# 1. What does one row look like?
var_summarized.iloc[0]
#important columns:GeneID, Type, ClinicalSignificance, PhenotypeList, OriginSimple, Chromosome

# 2. What are the unique clinical significance labels
var_summarized["ClinicalSignificance"].unique()

# 3. How many rows have "Benign" or "Pathogenic"?
print(var_summarized[var_summarized["ClinicalSignificance"] == "Benign"].sum())
print(var_summarized[var_summarized["ClinicalSignificance"] == "Pathogenic"].sum())

# 4. How big is the dataset? - number of rows vs columns

# 5. What information is available? - 5 columns that are important to the project

# 6. What are we predicting?

# 7. How clean is the data?

# 8. What types of variants exist?

# 9. How imbalanced is your target?

# 10. What does one single row look like?

# 11.  Where are the variants located?