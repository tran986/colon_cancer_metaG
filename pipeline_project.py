print("hello, I am still working on this project <3 ")
import urllib.request
import pandas as pd
import random
import os
import subprocess

#------------------------------------------------S1:DOWNLOAD RAW SEQ. SAMPLES FROM ENA: PRJEB7774
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
    fastq_ftp = seq_single[seq_single["sample_title"] == i]["fastq_ftp"]
    all_url=random.sample(fastq_ftp.tolist(), 10) #changes this if need more/less samples
    combined_url.extend(all_url)

# Directory to save files - Dowloads
out_dir = "ena_fastq"
os.makedirs(out_dir, exist_ok=True)

for url in combined_url:
    # Extract filename from URL
    filename = url.split("/")[-1]
    url = "https://" + url
    filepath = os.path.join(out_dir, filename)
    print(f"Downloading {url} -> {filepath}")
#   urllib.request.urlretrieve(url, filepath)

#make a metadata for downstream analysis later on:
metadata_samp = seq_info[seq_info["fastq_ftp"].isin(combined_url)][["run_accession","sample_title"]]
new_i=[]
for i in metadata_samp["sample_title"]:
    if i=="Stool sample from controls":
       new_i.append("Controls")
    else:
       if i=="Stool sample from carcinoma":
          new_i.append("Carcinoma")
       else:
          new_i.append("Advanced Adenoma")

metadata_samp["group_title"] = new_i
print(metadata_samp)

#------------------------------------------------S2: RUN QUALITY CONTROLS ON RAW SEQ:
#--------FASTQC:
qc_dir = "fastqc_dir"
print(os.listdir(out_dir))

#create qc_dir before running fastqc:
os.makedirs(qc_dir, exist_ok=True)

for filename in os.listdir(out_dir):
   filepath=os.path.join(out_dir, filename)
   subprocess.run(
   ["fastqc", filepath, "-o", qc_dir]
   )

#--------FASTP for trimming
trimmed_dir = "fastp_dir"
#create trimmed_dir before trimming:
os.makedirs(trimmed_dir, exist_ok=True)

for filename in os.listdir(out_dir):
   filepath_in=os.path.join(out_dir, filename)
   filepath_out=os.path.join(trimmed_dir, filename)
   print(f"trimming {filepath_in} to {filepath_out}")
   subprocess(
    "fastp",
    "-i", filepath_in,
    "-o", filepath_out,
    "--html", filepath_out + ".html"
   )
