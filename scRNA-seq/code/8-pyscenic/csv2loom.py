import os, sys 
os.getcwd() 
os.listdir(os.getcwd()) 
import loompy as lp; 
import numpy as np; 
import scanpy as sc; 
x=sc.read_csv("/home/huyifeng/pyscenic/Breg.all.csv"); ## 曾老师的代码这里是x=sc.read_csv("pbmc_3k.csv"); 
row_attrs = {"Gene": np.array(x.var_names),}; 
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("Breg.loom",x.X.transpose(),row_attrs,col_attrs);
