import csv
import pandas


transcript_list = []
protein1_list = []
protein2_list = []
protein3_list = []

for i in range(1, 151):
    transcript_list.append([i, "protein1", (0.15 * i)])
    transcript_list.append([i, "protein2", (0.10 * i)])
    transcript_list.append([i, "protein3", (0.06 * i)])
    if i == 150:
        protein1_list.append(0.15 * i)
        protein2_list.append(0.10 * i)
        protein3_list.append(0.06 * i)
for i in range(151, 301):
    transcript_list.append([i, "protein1", protein1_list[-1]])
    transcript_list.append([i, "protein2", (0.10 * i)])
    transcript_list.append([i, "protein3", (0.06 * i)])
# for i in range(226, 301):
#     transcript_list.append([i, "protein1", (-0.05 * i + 26.3)])
#     transcript_list.append([i, "protein2", protein2_list[-1]])
#     transcript_list.append([i, "protein3", (0.05 * i)])


file = open("../../../data/paper_data17.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
