from tabulate import tabulate
import os

name1 = ["3500_3503_","6995_3503_","6995_7031_"]
name2 = ["2","4","8","16"]
firstCol = ["3500:3503","6995:3503","6995:7031"]
head = ["SeqA:SeqB","Number of threads","Execution time(ms)","Speedup(compared to two threads)"]
table = [head]
for i in range(3):
    for j in range(4):
        temp = []
        temp.append(firstCol[i])
        temp.append(name2[j])
        path = "output/" + name1[i] + name2[j] + ".txt"
        with open(path,'r') as f:
            temp.append(int(f.readline()))
            temp.append(float(f.readline()))
        table.append(temp)
print(tabulate(table, headers='firstrow', tablefmt='fancy_grid', showindex=True))
os.system("pause")