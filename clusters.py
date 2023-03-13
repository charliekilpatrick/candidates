from astropy.table import Table
from astropy.io import fits
import numpy as np, os, sys, warnings, glob

def transpose(data):
    return(list(map(list, zip(*data))))

def import_cluster_file(filename):
        idx=0
        table_count=0
        filedata=[]
        photdata=[]
        canddata=[]
        phottabs={}
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('#'):
                    if table_count==0:
                        filehead=line.replace('#','').split()
                    elif table_count==1:
                        candhead=line.replace('#','').split()
                    else:
                        phothead=line.replace('#','').split()
                    if photdata:
                        phottabs[str(photdata[0][0].strip())]=photdata
                        photdata=[]
                    table_count += 1
                    continue
                if table_count==1:
                    filedata.append(line.split())
                elif table_count==2:
                    canddata.append(line.split())
                else:
                    photdata.append(line.split())
            if photdata:
                phottabs[str(photdata[0][0].strip())]=photdata
                photdata=[]

        filetable = Table(transpose(filedata), names=filehead)
        candtable = Table(transpose(canddata), names=candhead)
        for key in phottabs.keys():
            phottabs[key] = Table(transpose(phottabs[key]), names=phothead)
        return(filetable, candtable, phottabs)

def get_all_tables(table_dir):
    files = None
    candidates = None
    photometry = []
    for file in glob.glob(table_dir+'/*.clusters'):
        filetable, candtable, phottabs = import_cluster_file(file)
        if not files:
            files = Table(filetable)
        else:
            files = vstack([files, filetable])
        if not candidates:
            candidates = Table(candtable)
        else:
            candidates = vstack([candidates, candtable])
        photometry.append(phottabs)
    return(files, candidates, photometry)

files, candidates, photometry = get_all_tables('/data/LCO/Swope/workstch/gw190814_gw190814tmpl/1')
