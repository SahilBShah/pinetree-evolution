import csv
import pandas


transcript_list = []
protein1_list = []
protein2_list = []
protein3_list = []

for i in range(1, 301):
    transcript_list.append([i, "protein1", (0.12 * i)])
    transcript_list.append([i, "protein2", (0.11 * i)])
    transcript_list.append([i, "protein3", (0.10 * i)])
    transcript_list.append([i, "protein4", (0.09 * i)])
    transcript_list.append([i, "protein5", (0.08 * i)])
    transcript_list.append([i, "protein6", (0.07 * i)])
    transcript_list.append([i, "protein7", (0.06 * i)])
    transcript_list.append([i, "protein8", (0.05 * i)])
    transcript_list.append([i, "protein9", (0.04 * i)])
    transcript_list.append([i, "protein10", (0.03 * i)])
    # if i == 100:
    #     protein1_list.append(0.12 * i)
    #     protein2_list.append(0.07 * i)
    #     protein3_list.append(0.04 * i)
# for i in range(101, 251):
#     transcript_list.append([i, "protein1", protein1_list[-1]])
#     transcript_list.append([i, "protein2", protein2_list[-1]])
#     transcript_list.append([i, "protein3", (0.04 * i)])
# for i in range(251, 301):
#     transcript_list.append([i, "protein1", (-0.05 * i + 24.5)])
#     transcript_list.append([i, "protein2", protein2_list[-1]])
#     transcript_list.append([i, "protein3", (0.04 * i)])


file = open("../../../data/targets/paper_data145.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
