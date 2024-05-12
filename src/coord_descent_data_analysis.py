import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("finalproject.csv")
sea = sns.FacetGrid(df, col = "n_snvs")
sns_plot = sea.map(sns.scatterplot, "instance", "err")
sns_plot.figure.savefig("finalproject.png")