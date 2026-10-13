import json
import matplotlib.pyplot as plt

# Change this to the profile file to plot
fname = "profiles/profiles_Arthropoda_k6_w250_s20.json"

with open(fname) as f:
    prof = json.load(f)

# Sort by frequency
items = sorted(prof.items(), key=lambda x: -x[1])[30]

kmers = [x[0] for x in items]
vals  = [x[1] for x in items]

plt.figure(figsize=(10,4))
plt.bar(kmers, vals)
plt.xticks(rotation=90)
plt.ylabel("Frequency")
plt.title("Top 30 k-mers - Arthropoda (k=6)")
plt.tight_layout()
plt.show()
