import sys
import pysam
import pandas as pd

fin, fout = sys.argv[1:3]
df = pd.DataFrame(0, columns=['read_len', 'mapping_len'], index=range(0, 151))
samfile = pysam.AlignmentFile(fin, "rb")
align_len=pd.DataFrame(None, columns=["read_len", "mapping_len"])
for read in samfile:
    df.loc[read.query_length, 'read_len'] += 1
    df.loc[read.query_alignment_length, 'mapping_len'] += 1
samfile.close()
df.to_csv(fout)