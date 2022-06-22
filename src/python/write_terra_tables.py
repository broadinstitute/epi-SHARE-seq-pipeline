import csv
import argparse
import os 

class SubTableRow:
	def __init__(self, pkr, lib, r1, fq1, fq2, notes):
		self.pkr = pkr
		self.lib = lib
		self.r1 = r1
		self.fq1 = fq1
		self.fq2 = fq2
		self.notes = notes

class MainTableRow:
	def __init__(self, pkr, r1, rna_lib = None, rna_fq1 = None, rna_fq2 = None, atac_lib = None, atac_fq1 = None, atac_fq2 = None):
		self.pkr = pkr
		self.r1 = r1
		self.rna_lib = rna_lib if rna_lib is not None else []
		self.rna_fq1 = rna_fq1 if rna_fq1 is not None else []
		self.rna_fq2 = rna_fq2 if rna_fq2 is not None else []
		self.atac_lib = atac_lib if atac_lib is not None else []
		self.atac_fq1 = atac_fq1 if atac_fq1 is not None else []
		self.atac_fq2 = atac_fq2 if atac_fq2 is not None else []
		self.notes = notes

def make_sub_key(lib, r1, bcl):
	return lib + "-" + r1 + "-" + bcl

def make_main_key(pkr,r1):
	key = pkr + "-" + r1
	return key.replace(' ', '-')

def get_subset_names(fqs):
	basenames = [os.path.basename(x) for x in fqs]
	# Subset name is between Lane number and Read number
	return [x.split('_L')[1].split('_R')[0].split('_')[1] for x in basenames]

def update_sub_dict(key, pkr, lib, r1, fq1, fq2, typ, notes):
	if typ == 'ATAC':
		if key not in atac:
			atac[key] = SubTableRow(pkr, lib, r1, [fq1], [fq2], notes)
		else:
			atac[key].fq1.append(fq1)
			atac[key].fq2.append(fq2)
	if typ == 'RNA':
		if key not in rna:
			rna[key] = SubTableRow(pkr, lib, r1, [fq1], [fq2], notes)
		else:
			rna[key].fq1.append(fq1)
			rna[key].fq2.append(fq2)

def update_main_dict(key, pkr, r1, rna_lib = None, rna_fq1 = None, rna_fq2 = None, atac_lib = None, atac_fq1 = None, atac_fq2 = None):
	if rna_lib is None:
		if key not in main:
			main[key] = MainTableRow(pkr, [r1], atac_lib = [atac_lib], atac_fq1 = atac_fq1, atac_fq2 = atac_fq2)
		else:
			main[key].atac_lib.append(atac_lib)
			main[key].atac_fq1 += atac_fq1
			main[key].atac_fq2 += atac_fq2
	if atac_lib is None:
		if key not in main:
			main[key] = MainTableRow(pkr, [r1], rna_lib = [rna_lib], rna_fq1 = rna_fq1, rna_fq2 = rna_fq2)
		else:
			main[key].rna_lib.append(rna_lib)
			main[key].rna_fq1 += rna_fq1
			main[key].rna_fq2 += rna_fq2

def update_main_subsets(key, pkr, r1, rna_lib = None, rna_fq1 = None, rna_fq2 = None, atac_lib = None, atac_fq1 = None, atac_fq2 = None):
	if len(rna_lib) == 0:
		if key not in main:
			main[key] = MainTableRow(pkr, r1, atac_lib = atac_lib, atac_fq1 = atac_fq1, atac_fq2 = atac_fq2)
		else:
			main[key].r1 += r1
			main[key].atac_lib += atac_lib
			main[key].atac_fq1 += atac_fq1
			main[key].atac_fq2 += atac_fq2
	elif len(atac_lib) == 0:
		if key not in main:
			main[key] = MainTableRow(pkr, r1, rna_lib = rna_lib, rna_fq1 = rna_fq1, rna_fq2 = rna_fq2)
		else:
			main[key].r1 += r1
			main[key].rna_lib += rna_lib
			main[key].rna_fq1 += rna_fq1
			main[key].rna_fq2 += rna_fq2
	else:
		if key not in main:
			main[key] = MainTableRow(pkr, r1, rna_lib = rna_lib, rna_fq1 = rna_fq1, rna_fq2 = rna_fq2, atac_lib = atac_lib, atac_fq1 = atac_fq1, atac_fq2 = atac_fq2)
		else:
			main[key].r1 += r1
			main[key].atac_lib += atac_lib
			main[key].atac_fq1 += atac_fq1
			main[key].atac_fq2 += atac_fq2
			main[key].rna_lib += rna_lib
			main[key].rna_fq1 += rna_fq1
			main[key].rna_fq2 += rna_fq2

parser = argparse.ArgumentParser( description='Generate tables for Terra')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-n', '--name', type=str, required=True)
args = parser.parse_args()

atac = dict()
rna = dict()
main = dict()

in_fh = open(args.input)#args.input)
csv_reader = csv.DictReader(in_fh, delimiter='\t')

for record in csv_reader:
	pkr = record['PKR']
	lib = record['Library']
	fq1 = sorted(record['fastq_R1'].split(','))
	fq2 = sorted(record['fastq_R2'].split(','))
	r1 = get_subset_names(fq1)
	typ = record['Type']
	notes = record['Notes']
	for i in range(len(r1)):
		key = make_sub_key(lib, r1[i], args.name) #put args.name back in
		update_sub_dict(key, pkr, lib, r1[i], fq1[i], fq2[i], typ, notes)

with open('rna.tsv', 'wt') as outfile:
	tsv_writer = csv.writer(outfile, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
	tsv_writer.writerow(['entity:scRNA_libraries_id', 'BCL', 'PKR', 'Library', 'R1_Subset', 'fastq_R1', 'fastq_R2', 'Notes'])
	for key in rna.keys():
		pkr = rna[key].pkr
		lib = rna[key].lib
		r1 = rna[key].r1
		fq1 = list(rna[key].fq1)
		fq2 = list(rna[key].fq2)
		tsv_writer.writerow([key, args.name, pkr, lib, r1, '["' + '","'.join(fq1) + '"]', '["' + '","'.join(fq2) + '"]', rna[key].notes])
		main_key = make_main_key(pkr, r1)
		update_main_dict(main_key, pkr, r1, rna_lib=lib, rna_fq1=fq1, rna_fq2=fq2)

with open('atac.tsv', 'wt') as outfile:
	tsv_writer = csv.writer(outfile, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
	tsv_writer.writerow(['entity:scATAC_libraries_id', 'BCL', 'PKR', 'Library', 'R1_Subset', 'fastq_R1', 'fastq_R2', 'Notes'])
	for key in atac.keys():
		pkr = atac[key].pkr
		lib = atac[key].lib
		r1 = atac[key].r1
		fq1 = list(atac[key].fq1)
		fq2 = list(atac[key].fq2)
		tsv_writer.writerow([key, args.name, pkr, lib, r1, '["' + '","'.join(fq1) + '"]', '["' + '","'.join(fq2) + '"]', atac[key].notes])
		main_key = make_main_key(pkr, r1)
		update_main_dict(main_key, pkr, r1, atac_lib=lib, atac_fq1=fq1, atac_fq2=fq2)

# Merge all subsets together
for key in list(main.keys()):
	pkr = main[key].pkr
	r1 = list(main[key].r1)
	rna_lib = list(main[key].rna_lib)
	rna_fq1 = list(main[key].rna_fq1)
	rna_fq2 = list(main[key].rna_fq2)
	atac_lib = list(main[key].atac_lib)
	atac_fq1 = list(main[key].atac_fq1)
	atac_fq2 = list(main[key].atac_fq2)
	update_main_subsets(pkr.replace(' ', '_'), pkr, r1, rna_lib=rna_lib, rna_fq1=rna_fq1, rna_fq2=rna_fq2, atac_lib=atac_lib, atac_fq1=atac_fq1, atac_fq2=atac_fq2)

with open('run.tsv', 'wt') as outfile:
	tsv_writer = csv.writer(outfile, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
	tsv_writer.writerow(['entity:' + args.name + '_id', 'PKR', 'R1_Subset', 'ATAC_Lib', 'ATAC_fastq_R1', 'ATAC_fastq_R2', 'RNA_Lib', 'RNA_fastq_R1', 'RNA_fastq_R2'])
	for key in sorted(main.keys()):
		print(key)
		pkr = main[key].pkr
		r1 = '["' + '","'.join(main[key].r1) + '"]'
		atac_lib = '["' + '","'.join(set(main[key].atac_lib)) + '"]'
		atac_fq1 = '["' + '","'.join(main[key].atac_fq1) + '"]'
		atac_fq2 = '["' + '","'.join(main[key].atac_fq2) + '"]'
		rna_lib = '["' + '","'.join(set(main[key].rna_lib)) + '"]'
		rna_fq1 = '["' + '","'.join(main[key].rna_fq1) + '"]'
		rna_fq2 = '["' + '","'.join(main[key].rna_fq2) + '"]'
		tsv_writer.writerow([key.replace(" ", "_"), pkr, r1, atac_lib, atac_fq1, atac_fq2, rna_lib, rna_fq1, rna_fq2])
