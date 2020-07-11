import csv
import pandas


transcript_list = []
protein1_list = []
protein2_list = []
protein3_list = []

for i in range(1, 126):
    transcript_list.append([i, "protein1", (0.20 * i)])
    transcript_list.append([i, "protein2", (0.10 * i)])
    transcript_list.append([i, "protein3", (0.05 * i)])
    if i == 125:
        protein1_list.append(0.20 * i)
        protein2_list.append(0.10 * i)
        protein3_list.append(0.05 * i)
for i in range(126, 301):
    transcript_list.append([i, "protein1", protein1_list[-1]])
    transcript_list.append([i, "protein2", (0.10 * i)])
    transcript_list.append([i, "protein3", (0.05 * i)])
# for i in range(226, 301):
#     transcript_list.append([i, "protein1", (-0.05 * i + 26.3)])
#     transcript_list.append([i, "protein2", protein2_list[-1]])
#     transcript_list.append([i, "protein3", (0.05 * i)])


file = open("../../../data/targets/paper_data17.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
