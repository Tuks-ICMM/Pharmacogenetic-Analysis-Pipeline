import pandas as pd
import statsmodels.api as sm
import numpy as np

def pd_read_fwf (geneDict):
    
    dataFrames = dict()
    for geneName, geneFile in geneDict.items():
        listie = list()
        rf = open(geneFile)
        for row in rf.readlines():
            listie.append(row.split())
        rf.close()
        listie = pd.DataFrame(np.asarray(listie))
        listie.columns = listie.iloc[0]
        listie = listie[1:]
        dataFrames[str(geneName)] = listie
        dataFrames[str(geneName)]['SNP'] = dataFrames[str(geneName)]['SNP'].str.replace(r'^[0-9]{1,2}:[0-9]+[A-Z]+-[A-Z]+;(rs[0-9]+)$', r'\1', regex=True)
        dataFrames[str(geneName)] = dataFrames[str(geneName)][~dataFrames[str(geneName)].SNP.str.contains(r"<CN[0-9]+>")].astype({"MAC": "int32", "NCHROBS": "int32"})
        #dataFrames[str(geneName)] = dataFrames[str(geneName)].pivot(index="SNP", columns='CLST', values='MAC')
        dataFrames[str(geneName)] = dataFrames[str(geneName)].rename(columns={"MAC": "A1_OBS"}, errors="raise")[["SNP", "A1_OBS", "CLST"]]
        dataFrames[str(geneName)]['A2_OBS'] = dataFrames[str(geneName)].apply(lambda row: row.NCHROBS - row.A1_OBS, axis=1)
        dataFrames[str(geneName)] = pd.pivot_table(dataFrames[str(geneName)], index=['SNP', "A1"], columns=['CLST'], aggfunc='sum')
    return dataFrames