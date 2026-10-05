import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("final_dataset.csv")

# --- histogram ---

columns = [
    "Rg_md",
    "density_md",
    "end_to_end_md"
]


plt.figure(figsize=(12, 10))

for i, col in enumerate(columns, 1):
    plt.subplot(3, 2, i)
    
    sns.histplot(df[col], kde=True, bins=20, color="cornflowerblue")
    
    plt.title(col)

plt.tight_layout()
plt.savefig("hists.png")
plt.close()
