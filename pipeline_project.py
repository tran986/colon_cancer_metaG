print("hello, I am still working on this project <3 ")
import urllib.request
import pandas as pd
import random
import os

#------------------------------------------------S1:DOWNLOAD SAMPLES FROM ENA: PRJEB7774
#***Gut microbiome development along the colorectal adenoma-carcinoma sequence (Feng et al, 2015)

seq_info=pd.read_csv("~/Desktop/project_2/CRC_seq_info.txt", sep="\t")

#access sample ID from metadata:
print(seq_info.head(5))
print(seq_info.columns)

#only choose single-sq samples:
seq_single=seq_info[seq_info['submitted_ftp'].str.contains('clean.single')]
print(seq_single)

#Randomly choose 10 controls, 10 advanced adenoma, 10 carcinoma samples:

c = ["Stool sample from controls",
     "Stool sample from advanced adenoma",
     "Stool sample from carcinoma"]

combined_url=[]
for i in c:
    subset = seq_single[seq_single["sample_title"] == i]["fastq_ftp"]
    all_url=random.sample(subset.tolist(), 10)
    combined_url.extend(all_url)

print(combined_url)

# Directory to save files
out_dir = "~/Desktop/project_2/ena_fastq"
os.makedirs(out_dir, exist_ok=True)

for url in combined_url:
    # Extract filename from URL
    filename = url.split("/")[-1]
    filepath = os.path.join(out_dir, filename)
    print(f"Downloading {url} -> {filepath}")
    urllib.request.urlretrieve(url, filepath)







