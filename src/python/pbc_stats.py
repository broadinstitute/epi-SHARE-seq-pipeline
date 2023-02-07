import sys

# Author Kundaje lab
# https://github.com/kundajelab/ENCODE_scatac/blob/master/workflow/scripts/pbc_stats.py
# Input QNAME sorted


def calc_pbc(in_sam, out_path):
    total_pairs = 0
    distinct_pairs = -1
    one_read_pairs = 0
    two_read_pairs = 0

    current_pair = None
    current_count = 0

    for al in in_sam:
        fields = al.strip().split('\t')
        flag = int(fields[1])
        rname = fields[2]
        pos = int(fields[3])
        pnext = int(fields[7])

        if not (flag & 35 == 35):
            continue

        pair = (rname, pos, pnext)
        if pair == current_pair:
            total_pairs += 1
            current_count += 1
        else:
            total_pairs += current_count
            distinct_pairs += 1
            if current_count == 1:
                one_read_pairs += 1
            elif current_count == 2:
                two_read_pairs += 1

            current_pair = pair
            current_count = 1

    total_pairs += current_count
    distinct_pairs += 1
    if current_count == 1:
        one_read_pairs += 1
    elif current_count == 2:
        two_read_pairs += 1

    nrf = distinct_pairs / total_pairs
    pbc1 = one_read_pairs / distinct_pairs
    pbc2 = one_read_pairs / two_read_pairs

    stats_str = "\t".join(str(i) for i in [
        total_pairs,
        distinct_pairs,
        one_read_pairs,
        two_read_pairs,
        nrf,
        pbc1,
        pbc2
    ])
    descr_str = "\t".join([
        "TotalReadPairs",
        "DistinctReadPairs",
        "OneReadPair",
        "TwoReadPairs",
        "NRF=Distinct/Total",
        "PBC1=OnePair/Distinct",
        "PBC2=OnePair/TwoPair"
    ])
    with open(out_path, 'w') as f:
        f.write(f"{descr_str}\n{stats_str}\n")

if __name__ == "__main__":
    qc_path = sys.argv[1]
    calc_pbc(sys.stdin, qc_path)

