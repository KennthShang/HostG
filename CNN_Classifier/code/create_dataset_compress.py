import os
import numpy as np






if __name__ == '__main__':
    Load_path = "int_val/"
    name_list = os.listdir(Load_path)
    name_list.sort()

    #cnt = 0
    data = []

    for name in name_list:
        read = np.genfromtxt(Load_path+name, delimiter=',')
        if read.reshape(-1).shape[0] == 1998:
            read = read.reshape(1, 1998)
        
        data.append(read)
        np.savetxt("dataset/"+name, read, delimiter=",", fmt='%d')

