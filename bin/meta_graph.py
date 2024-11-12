import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys

input_string = sys.argv[1]
meta_id = sys.argv[2]

df = pd.read_csv(input_string, sep='\t', header=0)

my_plot = sns.catplot(y = 'Species', x = 'Abundance(%)', hue = 'Classifier', data = df, kind = 'bar', palette="pastel", edgecolor=".6", dodge=True, height=6, aspect=1.5).set(title=f"{meta_id}")

plt.yticks(fontsize=8)

ax = my_plot.facet_axis(0, 0)
for c in ax.containers:
    labels = [f'{(v.get_width())}' for v in c]
    ax.bar_label(c,label_type='edge', padding=2, fontsize=5.3)


plt.savefig(f'{meta_id}_classifier_comparison_plot.png', bbox_inches='tight')