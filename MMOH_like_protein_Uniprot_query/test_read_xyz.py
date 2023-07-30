"""
parse x,y,x file
"""

file_link = r"C:\Users\lisap\Downloads\GEO-pos-1.xyz"
colmarker = []
with open(file_link, 'r') as f:
    for line in f:
        print(line.split()[0])
