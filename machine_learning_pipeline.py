import pandas as pd
import numpy as np
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
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
#print(var_filter[["Name","ClinicalSignificance","ReferenceAlleleVCF", "AlternateAlleleVCF","PositionVCF", "Chromosome"]].head(5))

#what is number of benign and number of pathogenic after filtering:
var_filter.shape
(var_filter["ClinicalSignificance"] == "Benign").sum()

#2. Label encode your target column (ClinicalSignificance → 0/1)
var_label=var_filter["ClinicalSignificance"]
condition=[]
for i in range(len(var_label)):
   if var_label.iloc[i] == "Pathogenic":
       condition.append(1)
   else:
       condition.append(0)
var_filter["ClinicalSignificant_bin"]=condition

#3. One-Hot Encoding for ReferenceAlleleVCF and AlternateAlleleVCF:
ref_vcf=var_filter["ReferenceAlleleVCF"]
alter_vcf=var_filter["AlternateAlleleVCF"]

def aa_convert_func(ref_or_alter, ref_aa):
    res=[]
    for i in range(len(ref_or_alter)):
        if (ref_or_alter.iloc[i] == ref_aa):
            res.append(1)
        else:
            res.append(0)
    return res

base=["A", "T", "G", "C"]
aa_dict_ref = {}
aa_dict_alt = {}
for aa in (base):
    aa_dict_ref[f"ref_{aa}"]=aa_convert_func(ref_vcf, aa)
    aa_dict_alt[f"alter_{aa}"]=aa_convert_func(alter_vcf, aa)
    
encoded_vcf=pd.concat([pd.DataFrame(aa_dict_alt),
                       pd.DataFrame(aa_dict_ref)], axis=1)
encoded_vcf=encoded_vcf.fillna(0)
print(encoded_vcf.isnull().sum()) 

# 4: Put all X into 1 df, and Y into 1 df:
# Y - outcome var:
y_df=var_filter[["ClinicalSignificant_bin"]]
print(y_df.head(5))

# X - predictor var A: using VCF:
encoded_vcf.head(5)

# Sanity check: if x and y have same number of rows? any missing values in x? classes balanced?
#print(len(encoded_vcf) == len(y_df))
#print(encoded_vcf.isna().sum())
#print(encoded_vcf.isnull().sum())
#print(len(y_df[y_df["ClinicalSignificance"] == "Benign"]) == len(y_df[y_df["ClinicalSignificance"] != "Benign"]))

# ---------- STEP 2: Train/Test Split:
#1. split the dataset for training (80%) and test (20%)
x_train, x_test, y_train, y_test = train_test_split(
    encoded_vcf, #x df
    y_df,  #y outcome
    test_size = 0.2, #20% saved for testing
    random_state=123, #set.seed
    stratify=y_df #keeps benign/pathogenic ratio same in both splits
)

#print(x_train.shape[0])
#print(x_test.shape[0])

# ---------- STEP 3: Model Training and Evaluating:
#find a good number of tree? - loop through the loop of n possible of tree, compute score 
n_list = [10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 1000]
scores = []
for n in n_list: 
    model = RandomForestClassifier(
        n_estimators=n,
        random_state=123,
        class_weight="balanced"
    )
    score = cross_val_score(model, x_train, y_train, cv=5, scoring="roc_auc").mean() #calculate for score
    #print(f"score for {n}: {score:.4f}")
    scores.append(score)

score_all_est=pd.DataFrame({"n_tree":n_list,
            "score":scores})


plt.plot(score_all_est["n_tree"], score_all_est["score"])
plt.xlabel("number of trees")
plt.ylabel("score")
plt.title("trees vs scores")
#plt.show()

#fit model:
model = RandomForestClassifier(n_estimators=500, random_state=42, 
                                class_weight="balanced")
model.fit(x_train, y_train)

# 3. Predict
y_pred = model.predict(x_test)
y_prob = model.predict_proba(x_test)[:, 1]

print(classification_report(y_test, y_pred, 
      target_names=["Benign", "Pathogenic"]))
print("ROC-AUC:", roc_auc_score(y_test, y_prob))

y_prob = model.predict_proba(x_test)[:, 1]  
fpr, tpr, thresholds = roc_curve(y_test, y_prob)
roc_auc = auc(fpr, tpr)
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, 
         color="blue", 
         lw=2, 
         label=f"ROC Curve (AUC = {roc_auc:.2f})")

plt.plot([0, 1], [0, 1], 
         color="red", 
         lw=2, 
         linestyle="--", 
         label="Random Guessing (AUC = 0.50)")

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate (Recall)")
plt.title("ROC Curve — Variant Effect Predictor")
plt.legend(loc="lower right")
plt.grid(True)
plt.tight_layout()
plt.savefig("roc_curve.png", dpi=300)   # saves to your project folder
#plt.show()

# ---------- STEP 4: Model Training 2 + Add more features 
#take some more features to predict from var_filter:- Chromosome, 
print(var_filter.columns)
print(var_filter[["Chromosome","Type", "PositionVCF", "HGNC_ID"]])
print(var_filter.iloc[0])

var_filter["kabuki_pheno"] = np.where(var_filter["PhenotypeList"] == "Kabuki syndrome", 1, 0)
#Chromosome, PositionVCF, kabuki_pheno features added: - Put all X into 1 df, and Y into 1 df:
x_kabuki_df = var_filter[["kabuki_pheno", "Chromosome", "PositionVCF"]]

#clean up "Chromosome" column:
mask = x_kabuki_df["Chromosome"] != "na"
x_kabuki_df = x_kabuki_df[mask]
y_kabuki_df = y_df[mask] 
x_kabuki_df["Chromosome"]=x_kabuki_df["Chromosome"].replace({"X":'23'})
print(x_kabuki_df["Chromosome"].unique())

#0. make sure that x and y have the same number of rows
print(x_kabuki_df.shape)
print(y_kabuki_df.shape)

#1. split the dataset for training (80%) and test (20%)
x_train_2, x_test_2, y_train_2, y_test_2 = train_test_split(
    x_kabuki_df, #x df
    y_kabuki_df,  #y outcome
    test_size = 0.2, #20% saved for testing
    random_state=123, #set.seed
    stratify=y_kabuki_df #keeps benign/pathogenic ratio same in both splits
)
#2.fit model:
model.fit(x_train_2, y_train_2)

# 3. Predict
y_pred_2 = model.predict(x_test_2)
y_prob_2 = model.predict_proba(x_test_2)[:, 1]

print(classification_report(y_test_2, y_pred_2, 
      target_names=["Benign", "Pathogenic"]))
print("ROC-AUC:", roc_auc_score(y_test_2, y_prob_2))
"""
#4: extract model:
y_prob_2 = model.predict_proba(x_test_2)[:, 1]  
fpr, tpr, thresholds = roc_curve(y_test, y_prob_2)
roc_auc = auc(fpr, tpr)
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, 
         color="blue", 
         lw=2, 
         label=f"ROC Curve (AUC = {roc_auc:.2f})")

plt.plot([0, 1], [0, 1], 
         color="red", 
         lw=2, 
         linestyle="--", 
         label="Random Guessing (AUC = 0.50)")

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate (Recall)")
plt.title("ROC Curve — Variant Effect Predictor")
plt.legend(loc="lower right")
plt.grid(True)
plt.tight_layout()
plt.savefig("roc_curve.png", dpi=300)   # saves to your project folder
#plt.show()
"""
