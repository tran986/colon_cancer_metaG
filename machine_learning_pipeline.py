import pandas as pd
import numpy as np
# ---------- STEP 0: import ClinVar (NCBI dataset) and exploring the dataset: -------

# 1. Focusing on Benign and Pathogenic "only"
# 2. Filter out those located in germline and germline/somatic only
# 3. Select Type with SNV only:
"""
# ---------- STEP 1: Cleaning the data:
var_summarized=pd.read_csv("~/Desktop/colon_cancer_metaG/variant_summary.txt", sep = "\t", dtype={"Chromosome": str})
var_summarized=var_summarized.loc[
    (var_summarized["ClinicalSignificance"].isin(["Benign", "Pathogenic"])) &
    (var_summarized["OriginSimple"].isin(["germline","germline/somatic","not applicable"])) &
    (var_summarized["Type"] == "single nucleotide variant")
]

#drop columns that are NA
var_summarized=var_summarized.drop(columns=["ReferenceAllele", "AlternateAllele"])

# Subset the first 2k samples total (then combine)
var_filter=pd.concat([var_summarized[var_summarized["ClinicalSignificance"]=="Pathogenic"].head(1000),
                 var_summarized[var_summarized["ClinicalSignificance"]!="Pathogenic"].head(1000)], ignore_index = True)

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
   
#var_filter.to_csv("~/Desktop/colon_cancer_metaG/variant_subset.csv", index = False)
"""

var_filter=pd.read_csv("~/Desktop/colon_cancer_metaG/variant_subset.csv")

# 1. What does one row look like?
var_filter.iloc[0]

# 2. What are the unique clinical significance labels
var_filter["ClinicalSignificance"].unique()

# 3. How many rows have "Benign" or "Pathogenic"?
print(var_filter[var_filter["ClinicalSignificance"] == "Benign"].value_counts())
print(var_filter[var_filter["ClinicalSignificance"] == "Pathogenic"].value_counts())

# 4. How big is the dataset? - number of rows vs columns
var_filter.shape #(8907105 rows, 43 columns)

# 5. What information is available? - 5 columns that are important to the project
#important columns:GeneID, Type, ClinicalSignificance, PhenotypeList, OriginSimple, Chromosome
var_summarized_subset=var_filter[["GeneID", "SCVsForAggregateOncogenicityClassification", "ClinicalSignificance", "PhenotypeList", "OriginSimple", "PositionVCF"]]

# 6. What are we predicting?- predictor vars from the FULL dataset:
outcome_var=var_filter["ClinicalSignificance"]

# 7. How clean is the data?
#print(var_summarized.isna().values.any()) #check for NA - no NA values
var_filter.isnull().any #no NA values

# 8. What types of variants exist?
var_filter["Type"].unique()

# 9. How imbalanced is your target?
var_filter["Type"].value_counts() #snv is very dominant and complex is only 1

# 10.  Where are the variants located?
print(var_filter.columns)
var_filter["OriginSimple"].unique() #'germline', 'not applicable', 'germline/somatic', 'unknown'

# ---------- STEP 2: Feature Engineering ← k-mers, encode categories
# 1. Look at raw sequences:
#6 columns that are important ClinicalSignificance, ReferenceAlleleVCF, AlternateAlleleVCF, PositionVCF, Chromosome
print(var_filter[["Name","ClinicalSignificance","ReferenceAlleleVCF", "AlternateAlleleVCF","PositionVCF", "Chromosome"]].head(5))

#what is number of benign and number of pathogenic after filtering:
var_filter.shape
(var_filter["ClinicalSignificance"] == "Benign").sum()

# Label encode your target column (ClinicalSignificance → 0/1)
var_label=var_filter["ClinicalSignificance"]
condition=[]
for i in range(len(var_label)):
   if var_label.iloc[i] == "Pathogenic":
       condition.append(1)
   else:
       condition.append(0)
var_filter["ClinicalSignificant_bin"]=condition

#One-Hot Encoding for ReferenceAlleleVCF and AlternateAlleleVCF:
ref_vcf=var_filter["ReferenceAlleleVCF"]
alter_vcf=var_filter["AlternateAlleleVCF"]

def aa_convert_func(ref_or_alter, ref_aa):
    res=[]
    for i in range(len(ref_or_alter)):
        if (ref_or_alter[i] == ref_aa):
            res.append(1)
        else:
            res.append(0)
    return res

ref_list=["ref_A", "ref_T", "ref_G", "ref_C"]
aa_dict_ref = {}
aa_dict_alt = {}
for i in (ref_list):
    aa = i[4]
    aa_dict_ref[f"ref_{aa}"]=aa_convert_func(ref_vcf, aa)
    aa_dict_alt[f"alter_{aa}"]=aa_convert_func(alter_vcf, aa)
    
encoded_vcf=pd.concat([pd.DataFrame(aa_dict_alt),
                       pd.DataFrame(aa_dict_ref)])
encoded_vcf=encoded_vcf.fillna(0)




"""
Original:          One-Hot Encoded:
                   ref_A  ref_T  ref_C  ref_G 
"A"        →         1      0      0      0
"T"        →         0      1      0      0
"C"        →         0      0      1      0
"G"        →         0      0      0      1
Do the same for altA

Week 1:
 Step 1 — Label encode your target column (ClinicalSignificance → 0/1)
 Step 2 — One-hot encode ReferenceAllele and AlternateAllele
 Step 3 — Combine into one feature dataframe

Week 2:
 Step 4 — Write the get_kmers() function and test it on toy sequences
 Step 5 — Apply it to your real sequence column
 Step 6 — Use CountVectorizer to turn kmers into a matrix

"""
