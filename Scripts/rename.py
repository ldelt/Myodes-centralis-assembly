import os

for filename in os.listdir():
	if '~' in filename:
		basename, suffix = filename.rsplit('.fastq', maxsplit=1)
		new_filename = f'{basename}-{suffix[1:]}.fastq'
		os.rename(filename, new_filename)
