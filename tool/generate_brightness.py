from settings import *
import os


# This method is used to generate the nigh time brightness index.
# It expects wrld-rad.data.txt and v4.inv to be in the /tmp/input dir.
# The resulting file v4.inv will be placed in the /tmp/input dir as well.
def run():
    night_file = open(INPUT_DIR + 'wrld-rad.data.txt', 'r')
    inv_file = open(INPUT_DIR + 'v4.inv', 'r')
    new_inv = open(INPUT_DIR + 'v4_tmp.inv', 'w+')
    i_j_dict = {}
    for line in night_file:
        line = line.split()
        i = line[0]
        j = line[1]
        i_j = (i, j)
        i_j_dict[i_j] = line[2]

    for line in inv_file:
        inv_line = line.split()
        lon = inv_line[2]
        lat = inv_line[1]
        search_i = str(round((float(lon) + 180) * 120 + 1))
        search_j = str(round((21600 + .5 - (float(lat) + 90) * 120)))
        if int(search_j) >= 21600:
            search_j = '21600'
        if int(search_i) >= 43200:
            search_i = '1'
        try:
            new_inv.write(line.replace('\n', ' ' + i_j_dict[(search_i, search_j)]) + "     " + '\n')
        except:
            new_inv.write(line.replace('\n', ' ' + '0' + '     ' + '\n'))
    inv_file.close()
    new_inv.close()
    os.remove(INPUT_DIR + "v4.inv")
    os.rename(INPUT_DIR + "v4_tmp.inv", INPUT_DIR + "v4.inv")
