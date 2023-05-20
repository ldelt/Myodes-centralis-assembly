from Bio import SeqIO
import time

start = time.time() 

seen = set()
records = []

for record in SeqIO.parse("resulted_fail.fastq", "fastq"):  
    if record.seq not in seen:
        seen.add(record.seq)
        records.append(record)


#writing to a fasta file
SeqIO.write(records, "resulted_fail_wd.fastq", "fastq")
end = time.time()

print(f"Run time is {(end- start)/60}")
