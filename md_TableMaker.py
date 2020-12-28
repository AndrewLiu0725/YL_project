f = open("url.txt", "r+")
links = f.readlines()
f.close()
f = open("md.txt", "w+")
for link_index, link in enumerate(links):
    if link_index == 0:
        f.write("|"+link[:-1]+"|")
    elif link_index == 1:
        f.write(link[:-1]+"|\n")
        f.write("|:---:|:---:|\n")
    elif link_index % 2 == 0:
        f.write("|"+link[:-1]+"|")
    elif link_index % 2 == 1:
        f.write(link[:-1]+"|\n")
f.close()