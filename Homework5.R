# Reading and preprocessing the data
data = read.table("/home/raf/Survival\ analysis/Oscar.txt", header = TRUE, sep = "\t")
data['Age'] <- data['Final']-data['Birth']
