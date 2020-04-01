import csv
import pandas


transcript_list = []
proteinX_list = []
proteinY_list = []
proteinZ_list = []

for i in range(1, 251):
    transcript_list.append([i, "protein1", (0.16 * i)])
    if i == 150:
        proteinX_list = transcript_list[-1]
    transcript_list.append([i, "protein2", (0.12 * i)])
    if i == 150:
        proteinY_list = transcript_list[-1]
    transcript_list.append([i, "protein3", (0.08 * i)])
    if i == 150:
        proteinZ_list = transcript_list[-1]
# for i in range(151, 251):
#     transcript_list.append([i, "proteinX", (0.1 * i)])
#     transcript_list.append([i, "proteinY", proteinY_list[-1]])
#     transcript_list.append([i, "proteinZ", (-0.05 * i) + (0.15 * 200)])



file = open("../../data/grant_data/grant_data100.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
