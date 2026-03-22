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
var_summarized["OriginSimple"].unique() #'germline', 'not applicable', 'germline/somatic', 'unknown'

# ---------- STEP 1: Cleaning the data:
# 1. Focusing on Benign and Pathogenic "only"
# 2. Filter out those located in germline and germline/somatic only
# 3. Filter out those with SNP only:

var_filter=var_summarized.loc[
    (var_summarized["ClinicalSignificance"].isin(["Benign", "Pathogenic"])) &
    (var_summarized["OriginSimple"].isin(["germline","germline/somatic","not applicable"])) &
    (var_summarized["Type"] == "single nucleotide variant")
]
var_filter=var_filter.drop(columns=["ReferenceAllele", "AlternateAllele"])
#write a function the loop through all columns, print out their unique vars and if number of na:
unique_var=[]
na_num =[]

for col in var_filter.columns:
    unique=len(var_filter[col].unique())
    num= (var_filter[col] == "na").sum()
    unique_var.append(unique)
    na_num.append(num)

print(pd.DataFrame({"column_name":var_filter.columns,
                    "unique variables" : unique_var,
                    "number of na" : na_num})) #use var filter moving forward - its clean
   
# ---------- STEP 2: Feature Engineering ← k-mers, encode categories
# 1. Look at raw sequences:



"""
Week 1:
 Step 1 — Label encode your target column (ClinicalSignificance → 0/1)
 Step 2 — One-hot encode ReferenceAllele and AlternateAllele
 Step 3 — Combine into one feature dataframe

Week 2:
 Step 4 — Write the get_kmers() function and test it on toy sequences
 Step 5 — Apply it to your real sequence column
 Step 6 — Use CountVectorizer to turn kmers into a matrix

"""
