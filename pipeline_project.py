import urllib.request
import pandas as pd
import random
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage
from collections import defaultdict
import subprocess
import multiqc
import shutil
import gzip
import statsmodels.formula.api as smf
from glob import glob
from tabulate import tabulate

#------------------------------------------------S1:DOWNLOAD RAW SEQ. SAMPLES FROM ENA: PRJEB7774
#***Gut microbiome development along the colorectal adenoma-carcinoma sequence (Feng et al, 2015)
"""
seq_info = pd.read_csv("/Users/Nghitran/Desktop/colon_cancer_metaG/CRC_seq_info.txt", sep="\t")
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
    all_url=random.sample(fastq_ftp.tolist(), 5) #changes this if need more/less samples
    #combined_url.extend(all_url)

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
metadata_samp.to_csv("sample_info.csv", index=False)


#------------------------------------------------S2: RUN QUALITY CONTROLS ON RAW SEQ:
#--------FASTQC:

qc_dir = "fastqc_dir"
print(os.listdir(out_dir))

#create qc_dir before running fastqc:
os.makedirs(qc_dir, exist_ok=True)

for filename in os.listdir(out_dir):
   filepath=os.path.join(out_dir, filename)
   print(f"running QC on {filename} to {filepath}") #run fastqc
   subprocess.run(                             
   ["fastqc", filepath, "-o", qc_dir])
   
#--------MULTIQC:
multiqc_dir = "multiqc_dir"
os.makedirs(multiqc, exist_ok=True)
multiqc.run(qc_dir)            #run multiqc

#remove fastqc/multiqc dir after multiqc results are generated:
shutil.rmtree(qc_dir)
print(f"{qc_dir} has been removed to clear disk space")
shutil.rmtree("multiqc_data")
#print("multiqc_data has been remove to clear disk space")

#--------FASTP for trimming
trimmed_dir = "fastp_dir"
#create trimmed_dir before trimming:
os.makedirs(trimmed_dir, exist_ok=True)

for filename in os.listdir(out_dir):
   filepath_in=os.path.join(out_dir, filename)
   filepath_out=os.path.join(trimmed_dir, filename)
   print(f"trimming {filepath_in} to {filepath_out}")   #run fastp
   subprocess.run(["fastp",
    "-i", filepath_in,
    "-o", filepath_out,
    "--html", filepath_out + ".html"])

   
#------------------------------------------------S3 RUN METAPHLAN3 
# FOR TAXA PROFILING :
#splitting into chunks:-in order to run in local machine with low RAM:
tmp_dir="fastq_chunks"

os.makedirs(tmp_dir, exist_ok=True)

lines_per_chunk = 4_000_000  
nproc = 1

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
taxa_prof_dir = "metaphlan_dir"
os.makedirs(taxa_prof_dir, exist_ok=True)

#metaphlan version 3 (compatible with python 3.9)
#metaphlan --install --bowtie2db "/home/tran986/mp_database_new"

for filename in os.listdir(tmp_dir):
   if filename.endswith(".fastq.gz"):
   #if filename == "ERR710415.fastq.gz" : #test with 1 sample
      input_path=os.path.join(tmp_dir, filename)
      output_txt=filename.replace(".fastq.gz","_profile.txt")
      output_path=os.path.join(taxa_prof_dir,output_txt)
      print(f"profilling taxa in chunk {input_path} to {output_path}")
      subprocess.run(["metaphlan",input_path,
        "--input_type", "fastq", 
        "--nproc", "1",
        "--bowtie2db", "/Users/Nghitran/colon_cancer_metaG/metaphlan_database", #"--bowtie2db "/home/tran986/mp_database_new"
        "--index", "mpa_v31_CHOCOPhlAn_201901",
        "-o", output_path], check=True) #change dir to the db

#merge chunks for each sample:
merged_chunks_dir=os.path.join(taxa_prof_dir, "merge")
os.makedirs(merged_chunks_dir, exist_ok=True)

groups = defaultdict(list)
for f in os.listdir(taxa_prof_dir):
    if f.endswith("_profile.txt") and "_chunk_" in f:
        sample = f.split("_chunk_")[0]
        groups[sample].append(os.path.join(taxa_prof_dir, f))

for sample, files in groups.items():
    dfs = []
    for f in files:
        # Read data (ignore comments)
        df = pd.read_csv(f, sep="\t", index_col=0, comment="#")
        # Convert all values to numeric (force errors to NaN)
        df = df.apply(pd.to_numeric, errors='coerce')
        dfs.append(df)

    # Concatenate along columns (each chunk as a column)
    merged_df = pd.concat(dfs, axis=1)

    # Average across chunks
    merged_df = merged_df.mean(axis=1).to_frame(name="relative_abundance")

    # Add fake NCBI_tax_id column so we have 3 columns
    merged_df.insert(0, "NCBI_tax_id", 0)

    # Reset index to turn clade_name from index to column
    merged_df.reset_index(inplace=True)
    merged_df.rename(columns={"index": "clade_name"}, inplace=True)

    # Write back
    out_path = os.path.join(merged_chunks_dir, f"{sample}_profile.txt")
    print(f"Writing merged chunks to {out_path}")
    with open(out_path, "w") as out:
        out.write(f"#mpa_v31_CHOCOPhlAn_201901\n")
        out.write(f"#/opt/homebrew/Caskroom/miniforge/base/envs/metaphlan3/bin/metaphlan\n")
        out.write(f"#SampleID\tMetaphlan_Analysis\n")
        out.write(f"#clade_name\tNCBI_tax_id\trelative_abundance\n")
        merged_df.to_csv(out, sep="\t", index=False, header=False)

#merging all metaphlan output sample tbls:
merge_sample_path = os.path.join(taxa_prof_dir, "all_sample_merged_file.txt") #out

merge_all_samp=[]
for filename in os.listdir(merged_chunks_dir):
   if filename.endswith(".txt"):
      merge_samp=os.path.join(merged_chunks_dir, filename)
      merge_all_samp.append(merge_samp)
      print(f"Merging MetaPhlan all sample counts from {merge_samp}")
print(*merge_all_samp)

subprocess.run([
     "merge_metaphlan_tables.py", *merge_all_samp,
     "-o", merge_sample_path], check=True)

"""
#------------------------------------------------S4 READ INTO COUNT OUTPUTS:
#merge_file_path=os.path.join(merge_sample_dir, "all_sample_merged_file.txt") # Replace with actual file path
# Read the MetaPhlAn table
merge_sample_path=os.path.join("metaphlan_dir", "all_sample_merged_file.txt")

df_count_out = pd.read_csv(merge_sample_path, sep="\t", comment="#")
print(df_count_out.columns) #Index(['clade_name', 'NCBI_tax_id', 'ERR688431_profile', 'ERR688467_profile', 'ERR688483_profile']

#match the ERR samplename with condition:
profile=df_count_out.filter(regex='profile').columns.tolist()
accession=[]
for i in profile:
   part=i.split("_",1)
   accession.append(part[0])
print(accession)

seq_info = pd.read_csv("/Users/Nghitran/Desktop/colon_cancer_metaG/CRC_seq_info.txt", sep="\t")
seq_info_short=seq_info[["run_accession","sample_title"]]
seq_info_filter=seq_info_short[seq_info_short["run_accession"].isin(accession)] #CONDITIONS OF ALL THE SAMPLES USED
print(f"getting the sample's condition from metadata {seq_info_filter}")

#------------------------------------------------S5 APPLY NORMALIZATION (Compositionally add up to 1) AND CREATE HEATMAP:
profile_count=df_count_out.filter(regex='profile')
df_normalized=profile_count.div(profile_count.sum(axis=0), axis=1)
print(df_normalized)

#add taxa col in the front of profile count:
taxa_col = df_count_out.clade_name
df_normalized.insert(0, 'taxa', taxa_col)
print(f'printing normalized data: {df_normalized}')

#separate taxa col into k__, p__, c__,...
split_normalized_df= df_normalized['taxa'].str.split('|', expand=True)
print("\nSpliting data (intermediate DataFrame):")
split_normalized_df.columns = ['Domain', 'Phylum', 'Clade', 'Order', 'Family', 'Genus', 'Species']
df_normalized = pd.concat([df_normalized, split_normalized_df], axis=1)
df_normalized_fix=df_normalized.iloc[:, 1:len(df_normalized.columns)]
print(f"splited df: {df_normalized_fix}")

#look at the species level first (filter out None in Species level):
print("filtering Null species:")
df_species=df_normalized_fix[df_normalized_fix['Species'].notnull()]
df_species=df_species.loc[:, ['Species'] + df_species.columns[df_species.columns.str.contains('ERR')].tolist()]

#---make a heatmap (on normalized data)
t = [item + "_profile" for item in seq_info_filter.run_accession.tolist()] #make the run_accession column matched with "ERR.._profile"
seq_info_filter.loc[:, 'sample'] = t
cond_dict = dict(zip(seq_info_filter['sample'], seq_info_filter['sample_title'])) #convert 2 columns into a dict
sample_cond_gr=pd.Series(cond_dict).reindex(df_species.columns)
group_colors = sample_cond_gr.map({'Stool sample from controls': "skyblue", 
                                   'Stool sample from advanced adenoma': "salmon", 
                                   'Stool sample from carcinoma': "red"})

#sort
df_species_sorted = df_species[sample_cond_gr.sort_values().index]
group_colors_sorted = group_colors[sample_cond_gr.sort_values().index]

df_species_sorted = df_species_sorted.set_index('Species')

hm_normalized_species=sns.clustermap(
    df_species_sorted,
    col_colors=group_colors_sorted,
    cmap="viridis",
    vmin=0.003,
    vmax=0.05,
    row_cluster=True,
    col_cluster=False,
    figsize=(7,6)
)

handles = [
    mpatches.Patch(color='skyblue', label='controls'),
    mpatches.Patch(color='salmon', label='advanced adenoma'),
    mpatches.Patch(color='red', label='carcinoma')
]
hm_normalized_species.ax_col_dendrogram.legend(
    handles=handles, 
    loc='center', 
    ncol=3, 
    bbox_to_anchor=(1.1, 0.72)
)
plt.show()

#------------------------------------------------APPLY LINEAR MODEL AND SIGNIFICANCE TEST(ANOVA) + VISUALIZATION
print(seq_info_filter)
print(df_species)

#pivot_longer the df_species: pivot_longer(col=!Species, names_to=sample, values_to=count)
df_species_longer=df_species.melt(id_vars='Species', 
                                  var_name='sample', 
                                  value_name='count')

#merge with condition:
df_species_merge=pd.merge(seq_info_filter, df_species_longer, on='sample', how='inner')
print(df_species_merge)

#fit linear model: --nolog:
model = smf.ols("count ~ sample_title + Species", data=df_species_merge).fit()
print(model.summary())
summary_df = pd.DataFrame({
    'term': model.params.index,
    'estimate': model.params.values,
    'p_value': model.pvalues.values,
    'std_err': model.bse.values,
    't_value': model.tvalues.values
})
#fit linear model -- log of relative abundance: result yields more hits due to compositionality problem addressed
df_species_merge["log_count"]=np.log1p(df_species_merge["count"])
model_log=smf.ols("log_count ~ sample_title + Species", data=df_species_merge).fit()
print(model_log.summary())

summary_df_log = pd.DataFrame({
    'term': model_log.params.index,
    'estimate': model_log.params.values,
    'p_value': model_log.pvalues.values,
    'std_err': model_log.bse.values,
    't_value': model_log.tvalues.values
})
print(summary_df)
print(summary_df_log)

#filter significant species (alpha level 0.05)
significant_tax_log=summary_df_log[summary_df_log.p_value <= 0.05]
print(significant_tax_log)
t=[]
for i in significant_tax_log.term:
   part=i.split(".")[1].rstrip("]")
   t.append(part)
print(t)

significant_tax_log["species"]=t
#export significant result:
print(tabulate(significant_tax_log[["species","estimate","p_value","std_err","t_value"]].reset_index(drop=True), headers='keys', tablefmt='psql'))

#--make a graph showing changes for these species 
print((significant_tax_log.species)) #17

bar_col = {'Stool sample from controls': 'blue', 
           'Stool sample from advanced adenoma': 'orange', 
           'Stool sample from carcinoma': 'red'}

for s in significant_tax_log["species"]:
   bar=sns.barplot(x='sample_title', 
            y='count', 
            data=df_species_merge[df_species_merge["Species"]==s],
            estimator="mean",
            errorbar="sd",
            capsize=0.2,
            palette=bar_col,
            order=['Stool sample from controls','Stool sample from advanced adenoma','Stool sample from carcinoma',])
   bar.set_xticklabels(['Control', 'Advanced Adenoma', 'Carcinoma'])
   plt.ylim(0, 0.03)
   plt.xlabel("Sample Group")
   plt.title(f"Normalized Count Across 3 Groups in {s}")
   plt.ylabel("Count")
   plt.show()

#------------------------------------------------BETA DIVERSITY + FEATURE REDUCTION (PCoA/PCA) + VISUALIZATION


#------------------------------------------------FUNCTIONAL ANALYSIS/NETWORKING:
