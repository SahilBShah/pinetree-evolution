import csv
import pandas


transcript_list = []
protein1_list = []
protein2_list = []
protein3_list = []

for i in range(1, 151):
    transcript_list.append([i, "protein1", (0.16 * i)])
    if i == 150:
        protein1_list = transcript_list[-1]
    transcript_list.append([i, "protein2", (0.12 * i)])
    if i == 150:
        protein2_list = transcript_list[-1]
    transcript_list.append([i, "protein3", (0.08 * i)])
    if i == 150:
        protein3_list = transcript_list[-1]
for i in range(151, 251):
    transcript_list.append([i, "protein1", protein1_list[-1]])
    transcript_list.append([i, "protein2", protein2_list[-1]])
    transcript_list.append([i, "protein3", protein3_list[-1]])
# for i in range(151, 251):
#     transcript_list.append([i, "proteinX", (0.1 * i)])
#     transcript_list.append([i, "proteinY", proteinY_list[-1]])
#     transcript_list.append([i, "proteinZ", (-0.05 * i) + (0.15 * 200)])



file = open("../../data/grant_data14.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
