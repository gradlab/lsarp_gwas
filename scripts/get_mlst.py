#!/usr/bin/env python3

batch_dict = {}
mlst_dict = {}
contaminated_samples = []

with open("data/contaminated_samples.txt", "r") as infile:
    for i,line in enumerate(infile):
        if i > 0:
            line = line.strip()
            contaminated_samples.append(line)

with open("data/sequenced_isolates.txt", "r") as infile:
    for i,line in enumerate(infile):
        if i > 0:
            line = line.strip().split()
            if line[1] not in contaminated_samples:
                batch_dict[line[1]] = line[0]

for samp in batch_dict:
    if batch_dict[samp] == "batch_2":
        infile_name = f"../../pipeline/Staphylococcus_aureus/batch_2/results/{samp}/MLST/ariba_mlst/mlst_report.tsv"
    else:
        infile_name = f"../../pipeline/Staphylococcus_aureus/{batch_dict[samp]}/final_results/{samp}/MLST/mlst_report.tsv"
    with open(infile_name, "r") as infile:
        for i,line in enumerate(infile):
            if i > 0:
                mlst_dict[samp] = line.strip().split()[0]
        
with open("data/population_structure/sample_mlst.txt", "w") as outfile:
    outfile.write("Sample\tST\n")
    for sample in mlst_dict:
        outfile.write(f"{sample}\t{mlst_dict[sample]}\n")
