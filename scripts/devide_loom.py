import loompy
import random
import numpy as np
import pandas as pd
import re

ds = loompy.connect('Raw/adult_human_20221007.loom')
annotation = pd.read_excel("Raw/cluster_annotation.xlsx")
map_dict = annotation.set_index("Cluster ID")["Supercluster"].to_dict()
cluster_values = ds.ca["Clusters"]
cellType_values = [map_dict.get(i, i) for i in cluster_values]
ds.ca["cellType"] = cellType_values

# Get the unique values of the 'ROIGroupCoarse' attribute
tissues = np.unique(ds.ca.ROIGroupCoarse)

for tissue in tissues:
    t = ds.ca.ROIGroupCoarse == tissue 
    # get the indices where ROIGroupCoarse == tissue
    indices = [i for i, x in enumerate(t) if x]
    # subsample indices 
    indices_sub = random.sample(indices, 5000)
    labels = ds.ca.cellType[np.sort(indices_sub)]
    clusters = ds.ca.Clusters[np.sort(indices_sub)]
    donor = ds.ca.CellID[np.sort(indices_sub)]
    cellID = ds.ca.Donor[np.sort(indices_sub)]
    counts_sub = ds[:, np.sort(indices_sub)]
    row_attrs = { "Gene": [re.sub(r'\.\d+$', '', item) for item in ds.ra.Accession]}
    col_attrs = {"CellID": cellID ,"cellType": labels, 'sampleID':donor, 'clusters': clusters }

    tissue = re.sub(" ", "", tissue)
    loompy.create('Split/' + tissue +'_5000.loom', counts_sub, row_attrs, col_attrs)

ds.close()