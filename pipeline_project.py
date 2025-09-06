print("hello, I am still working on this project <3 ")
import urllib.request
import pandas as pd
import random
import os
from collections import defaultdict
import subprocess
#import multiqc
import shutil
import gzip
from glob import glob

#------------------------------------------------S1:DOWNLOAD RAW SEQ. SAMPLES FROM ENA: PRJEB7774
#***Gut microbiome development along the colorectal adenoma-carcinoma sequence (Feng et al, 2015)

seq_info = pd.read_csv("/Users/Nghitran/colon_cancer_metaG/CRC_seq_info.txt", sep="\t")
#seq_info = pd.read_csv("/mnt/c/Users/Nia Tran/colon_cancer_metaG/CRC_seq_info.txt", sep="\t")

#access sample ID from metadata:
print(seq_info.head(5))
print(seq_info.columns)

#only choose single-sq samples:
seq_single=seq_info[seq_info['submitted_ftp'].str.contains('clean.single')]
print(seq_single)

#Randomly choose 1 controls, 1 advanced adenoma, 5 carcinoma samples:

c = ["Stool sample from controls",
     "Stool sample from advanced adenoma",
     "Stool sample from carcinoma"]

combined_url=[]
for i in c:
    fastq_ftp = seq_single[seq_single["sample_title"] == i]["fastq_ftp"]
    all_url=random.sample(fastq_ftp.tolist(), 1) #changes this if need more/less samples
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
    #urllib.request.urlretrieve(url, filepath)    #download fastq seq

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
   print(f"running QC on {filename} to {filepath}") #run fastqc
   """
   subprocess.run(                             
   ["fastqc", filepath, "-o", qc_dir]
   )
   """

#--------MULTIQC:
multiqc_dir = "multiqc_dir"
#multiqc.run(qc_dir, "-o", multiqc_dir)            #run multiqc

#remove fastqc/multiqc dir after multiqc results are generated:
#shutil.rmtree(qc_dir)
print(f"{qc_dir} has been removed to clear disk space")
#shutil.rmtree(multiqc_dir)
print(f"{multiqc_dir} has been remove to clear disk space")


#--------FASTP for trimming
trimmed_dir = "fastp_dir"
#create trimmed_dir before trimming:
os.makedirs(trimmed_dir, exist_ok=True)

for filename in os.listdir(out_dir):
   filepath_in=os.path.join(out_dir, filename)
   filepath_out=os.path.join(trimmed_dir, filename)
   print(f"trimming {filepath_in} to {filepath_out}")   #run fastp
   """
    subprocess.run([
    "fastp",
    "-i", filepath_in,
    "-o", filepath_out,
    "--html", filepath_out + ".html"])
   """


#------------------------------------------------S3 RUN METAPHLAN 4 FOR TAXA PROFILING :
#splitting into chunks:-in order to run in local machine with low RAM:
tmp_dir="fastq_chunks"

os.makedirs(tmp_dir, exist_ok=True)

lines_per_chunk = 4_000_000  
nproc = 1
"""
print("Splitting FASTQ into smaller chunks...")
for filename in os.listdir(trimmed_dir):
    if filename.endswith(".fastq.gz"):
        filepath = os.path.join(trimmed_dir, filename)
        sample_name=filename.split(".")[0]
        with gzip.open(filepath, "rt") as f: #open fastq file in txt mode
            chunk_idx = 0
            lines_buffer = []
            for i, line in enumerate(f, 1):
                lines_buffer.append(line)
                if i % lines_per_chunk == 0:
                    chunk_name = os.path.join(tmp_dir, f"{sample_name}_chunk_{chunk_idx}.fastq.gz")
                    with gzip.open(chunk_name, "wt") as cf:
                        cf.writelines(lines_buffer)
                        lines_buffer = []
                        chunk_idx += 1
    # Write remaining lines
            if lines_buffer:
                chunk_name = os.path.join(tmp_dir, f"{sample_name}_chunk_{chunk_idx}.fastq.gz")
                with open(chunk_name, "wt") as cf:
                    cf.writelines(lines_buffer)
                        
        print(f"Created {chunk_idx + 1} chunks for {sample_name}.")

print("Make sure they all gunzipped (insanity check):")

for filename in os.listdir(tmp_dir):
    filepath = os.path.join(tmp_dir, filename)
    # Skip directories
    if os.path.isdir(filepath):
        continue

    # Check if file is already gzipped
    with open(filepath, 'rb') as f:
        magic_number = f.read(2)

    if magic_number != b'\x1f\x8b':  # not gzipped
        if filepath.endswith(".fastq.gz"):
            new_path = filepath[:-3]  # remove .gz extension
            os.rename(filepath, new_path)       # rename all to unzipped
            filepath_gz=str(new_path) + ".gz"
            print(new_path)
            print(filepath_gz)
            with open(new_path, 'rb') as f_in, gzip.open(filepath_gz, 'wb') as f_out:
               f_out.writelines(f_in)
            print(f"Gzipped {filename}")
    else:
        print(f"Already gzipped: {filename}")


#run "rm *.fastq" if there are replicates:
for filename in os.listdir(tmp_dir):
    if filename.endswith(".fastq"):
        filepath=os.path.join(tmp_dir, filename)
        os.remove(filepath)
        print(f"Delete {filename} because of unzipp")

#run Metaphlan3:
"""
taxa_prof_dir = "metaphlan_dir"
os.makedirs(taxa_prof_dir, exist_ok=True)
"""
#metaphlan version 3 (compatible with python 3.9)
#metaphlan --install --bowtie2db "/home/tran986/mp_database_new"

for filename in os.listdir(tmp_dir):
   if filename.endswith(".fastq.gz"):
   #if filename == "ERR710415.fastq.gz" : #test with 1 sample
      input_path=os.path.join(tmp_dir, filename)
      output_txt=filename.replace(".fastq.gz","_profile.txt")
      output_path=os.path.join(taxa_prof_dir,output_txt)
      print(f"profilling taxa in chunk {input_path} to {output_path}")
      subprocess.run(["metaphlan",
                      input_path,
                      "--input_type", "fastq", 
                      "--nproc", "1",
                      "--bowtie2db", "/Users/Nghitran/colon_cancer_metaG/metaphlan_database", #"--bowtie2db "/home/tran986/mp_database_new"
                      "--index", "mpa_v31_CHOCOPhlAn_201901",
                      "-o", output_path  
                      ], check=True) #change dir to the db
      
"""
#merge chunks for each sample:
taxa_prac="metaphlan_dir_cp"
merged_chunks_dir="metaphlan_dir_cp/merge"
os.makedirs(merged_chunks_dir, exist_ok=True)
"""
groups = defaultdict(list)
for f in os.listdir(taxa_prac):
    if f.endswith("_profile.txt") and "_chunk_" in f:
        sample = f.split("_chunk_")[0]
        groups[sample].append(os.path.join(taxa_prac, f))

for sample, files in groups.items():
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep="\t", index_col=0, comment="#")
        dfs.append(df)

    # Concatenate along columns (each chunk as a column)
    merged_df = pd.concat(dfs, axis=1)

    # Average across chunks (to keep values in relative abundance scale)
    #merged_df = merged_df.mean(axis=1).to_frame(name=sample)

     # Write back with MetaPhlAn headers
    print("Writing sum of merged chunks to {merged_chunk_dir}")
    out_path = os.path.join(merged_chunks_dir, f"{sample}_profile.txt")
    with open(out_path, "w") as out:
        out.write(f"#SampleID\tMetaphlan_Analysis\n")
        out.write(f"#clade_name\tNCBI_tax_id\trelative_abundance\n")
        merged_df.to_csv(out, sep="\t")
        
"""
#merging all metaphlan output sample tbls:
merge_sample_dir = os.path.join(taxa_prac, "all_sample_merged_file.txt") #out
merged_chunks_dir=os.path.join(taxa_prac, "merge")
merge_all_samp=[]
for filename in os.listdir(merged_chunks_dir):
   if filename.endswith(".txt"):
      merge_samp=os.path.join(merged_chunks_dir, filename)
      merge_all_samp.append(merge_samp)
      print(f"Merging MetaPhlan all sample counts from {merge_samp}")
print(*merge_all_samp)

subprocess.run([
     "merge_metaphlan_tables.py", *merge_all_samp,
     "-o", merge_sample_dir], check=True)

