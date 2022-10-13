#!/usr/bin/env python
# coding: utf-8

# In[6]:


import subprocess
import sys
import os
import glob
if __name__ == '__main__':
    files = glob.glob("../Data/Diff_UV_Probes/*.tsv")
    for file in files:
        data_name = file.split('/')[-1][:-4]
        output_cluster = f"../Data/Diff_UV_Probes/{data_name}_cluster.tsv"
        output_matrix = f"../Data/Diff_UV_Probes/{data_name}_distance_matrix.npz"
        command_list = ["python",
                        "affinityProp.py",
                        "-i",
                        file,
                        "-oC",
                        output_cluster,
                        "-oM",
                        output_matrix,
                       "-c",
                       '5',
                       "-r",
                       '0']
        affinityProp = subprocess.Popen(command_list)
        affinityProp.communicate()


# In[ ]:




