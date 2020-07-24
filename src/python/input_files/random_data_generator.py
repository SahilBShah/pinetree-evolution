import csv
import pandas


transcript_list = []
protein1_list = []
protein2_list = []
protein3_list = []

for i in range(1, 101):
    transcript_list.append([i, "protein1", (0.04 * i)])
    transcript_list.append([i, "protein2", (0.07 * i)])
    transcript_list.append([i, "protein3", (0.12 * i)])
    if i == 100:
        protein1_list.append(0.04 * i)
        protein2_list.append(0.07 * i)
        protein3_list.append(0.12 * i)
for i in range(101, 251):
    transcript_list.append([i, "protein1", (0.04 * i)])
    transcript_list.append([i, "protein2", protein2_list[-1]])
    transcript_list.append([i, "protein3", protein3_list[-1]])
for i in range(251, 301):
    transcript_list.append([i, "protein1", (0.04 * i)])
    transcript_list.append([i, "protein2", protein2_list[-1]])
    transcript_list.append([i, "protein3", (-0.05 * i + 24.5)])


file = open("../../../data/targets/paper_data139.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
