import gzip
import sys

REV_COMP = str.maketrans("ATGC", "TACG")
def reverse_complement(seq):
    return str.translate(seq, REV_COMP)[::-1]

def get_open_fn(path):
    with open(path, "rb") as f:
        is_gzipped = (f.read(2) == b'\x1f\x8b')
    return gzip.open if is_gzipped else open

def read_barcodes(path):
    open_fn = get_open_fn(path)
    with open_fn(path, 'rt') as file:
        bc = [b.strip() for b in file]
    bcrc = [reverse_complement(b) for b in bc]
    return set(bc), set(bcrc)

def bc_detect(fastq, whitelist, out, qc, offset, num_reads=100000, thresh=0.45):
    bc, bcrc = read_barcodes(whitelist)

    bc_match = 0
    bcrc_match = 0
    num_lines = num_reads * 4
    with gzip.open(fastq, 'rt') as f:
        for lnum, line in enumerate(f):
            if lnum >= num_lines:
                break
            if lnum % 4 != 1:
                continue
            seq = line.strip()[offset:]
            if seq in bc:
                bc_match += 1
            if seq in bcrc:
                bcrc_match += 1

    bc_match_prop = bc_match / num_reads
    bcrc_match_prop = bcrc_match / num_reads
    valid = (bc_match_prop >= thresh) or (bcrc_match_prop >= thresh)
    fc_chosen = (bc_match_prop >= bcrc_match_prop)

    with open(qc, 'w') as f:
        f.write(f"Direct match proportion: {bc_match_prop}\n")
        f.write(f"Reverse-complement match proportion: {bcrc_match_prop}\n")
        f.write(f"Reverse-complement chosen: {not fc_chosen}\n")

    if not valid:
        raise ValueError(f"Insufficient barcode match rate: {bc_match_prop}, {bcrc_match_prop}")
    with open(out, 'w') as f:
        if fc_chosen:
            f.write(f"{0}\n")
        else:
            f.write(f"{1}\n")

try:
    fastq = sys.argv[1]
    modality = sys.argv[2]
    whitelist = sys.argv[3]

    qc = sys.argv[4]
    out = sys.argv[5]
    thres = sys.argv[6]

    if modality == "10x":
        offset = 0
        bc_detect(fastq, whitelist, out, qc, offset, 100000, float(thres))
    elif modality == "10x_multiome":
        offset = 8
        bc_detect(fastq, whitelist, out, qc, offset, 100000, float(thres))

except NameError:
    pass
