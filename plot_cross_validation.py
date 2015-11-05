import glob, json, collections
import pandas as pd

F_RESULTS = glob.glob("clf/*.json")
results = collections.defaultdict(list)

KWIN = [0,1,2,]

data = []
for f in F_RESULTS:
    with open(f) as FIN:
        js = json.load(FIN)
    js.pop("train_pdb")
    js.pop("test_pdb")
    data.append(js)

df = pd.DataFrame(data)
df.sort(columns=["kernel_window",],inplace=True)

import matplotlib.pyplot as plt
import seaborn as sns
fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
ax0,ax1 = axes.ravel()

text = "{}-fold cross-validation {} residues off diagonal"
text = text.format(
    df["k_fold_n"][0],
    df["diag_window"][0],
)
            
sns.set_style("whitegrid")
sns.barplot(data=df, x="kernel_window",y="test_acc",ax=ax0, hue="ratio_TP_to_FP")
sns.barplot(data=df, x="kernel_window",y="ROC_AUC",ax=ax1, hue="ratio_TP_to_FP")
ax0.set_ylim(0.6, 1.0)
ax0.set_title(text)
sns.plt.show()
