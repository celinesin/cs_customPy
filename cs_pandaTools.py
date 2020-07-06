import numpy as np
import pandas as pd
from scipy.spatial import distance

# calculates jaccard indices for on columns of dataframe and returns dataframe of values
# amended from pandas.corr() from frame.py
def ji(exprdata, min_periods=1):

    numeric_df = exprdata._get_numeric_data()
    cols = numeric_df.columns
    idx = cols.copy()
    mat = numeric_df.values.T

    if min_periods is None:
        min_periods = 1

    K = len(cols)
    jimat = np.empty((K, K), dtype=float)
    mask = np.isfinite(mat)

    for i, ac in enumerate(mat):
        for j, bc in enumerate(mat):
            if i > j:
                continue

            valid = mask[i] & mask[j]
            if valid.sum() < min_periods:
                c = np.nan
            elif i == j:
                c = 1.0
            elif not valid.all():
                c = distance.jaccard(ac[valid], bc[valid])
            else:
                c = distance.jaccard(ac, bc)
            jimat[i, j] = c
            jimat[j, i] = c

    return pd.DataFrame(jimat, index=idx, columns=cols)
