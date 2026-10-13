import matplotlib.pyplot as plt

# Read lengths
with open("lengths.txt") as f:
    lengths = [int(x.strip()) for x in f if x.strip()]

plt.hist(lengths, bins=80, color='steelblue')
plt.xlabel("Length (bp)")
plt.ylabel("Frequency")
plt.title("Sequence length distribution")
plt.tight_layout()
plt.show()
