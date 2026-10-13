import pandas as pd
import re

tax = pd.read_csv("data/raw/taxon.taxon", sep="\t", header=None, names=["seqid", "taxonomy"])

def get_phylum(s):
    m = re.search(r"p__([A-Za-z0-9]+)", s)
    return m.group(1) if m else None

tax["phylum"] = tax["taxonomy"].apply(get_phylum)

# list unique phyla

print("Unique phyla:") 
print(tax["phylum"].unique())

# count sequences per phylum
print("\nCount per phylum:")
print(tax["phylum"].value_counts())
