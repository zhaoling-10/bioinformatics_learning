#!/usr/bin/env python3
# paf_to_chain.py  ← save as this file name
# method: python3 paf_to_chain.py input.paf output.chain

import sys

def paf_to_chain(paf_file, chain_file):
    chain_id = 0
    with open(paf_file) as fin, open(chain_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            f = line.strip().split('\t')
            if len(f) < 12:
                continue

            # PAF fields
            qname  = f[0]
            qlen   = int(f[1])
            qstart = int(f[2])
            qend   = int(f[3])
            strand = f[4]
            tname  = f[5]
            tlen   = int(f[6])
            tstart = int(f[7])
            tend   = int(f[8])
            nmatch = int(f[9])
            score  = nmatch

            # Extract cs or cg string (for building gap information)
            cs_str = None
            for field in f[12:]:
                if field.startswith('cg:Z:'):
                    cs_str = field[5:]
                    break

            if cs_str is None:
                # No cigar, output simple single block chain
                chain_id += 1
                size = qend - qstart
                if strand == '+':
                    fout.write(f"chain {score} {tname} {tlen} + "
                               f"{tstart} {tend} "
                               f"{qname} {qlen} + "
                               f"{qstart} {qend} {chain_id}\n")
                else:
                    fout.write(f"chain {score} {tname} {tlen} + "
                               f"{tstart} {tend} "
                               f"{qname} {qlen} - "
                               f"{qlen-qend} {qlen-qstart} {chain_id}\n")
                fout.write(f"{size}\n\n")
                continue

            # Parse CIGAR to generate chain format
            import re
            ops = re.findall(r'(\d+)([MIDNSHP=X])', cs_str)
            
            chain_id += 1
            if strand == '+':
                fout.write(f"chain {score} {tname} {tlen} + "
                           f"{tstart} {tend} "
                           f"{qname} {qlen} + "
                           f"{qstart} {qend} {chain_id}\n")
            else:
                fout.write(f"chain {score} {tname} {tlen} + "
                           f"{tstart} {tend} "
                           f"{qname} {qlen} - "
                           f"{qlen-qend} {qlen-qstart} {chain_id}\n")

            blocks = []
            cur_match = 0
            dt = 0  # gap in target
            dq = 0  # gap in query

            for length, op in ops:
                length = int(length)
                if op in ('M', '=', 'X'):
                    if dt > 0 or dq > 0:
                        if cur_match > 0:
                            blocks.append((cur_match, dt, dq))
                        cur_match = length
                        dt = 0
                        dq = 0
                    else:
                        cur_match += length
                elif op == 'D':
                    dt += length
                elif op in ('I', 'S'):
                    dq += length
                elif op == 'N':
                    dt += length

            if cur_match > 0:
                blocks.append((cur_match, 0, 0))

            for i, (size, dt, dq) in enumerate(blocks):
                if i < len(blocks) - 1:
                    fout.write(f"{size}\t{dt}\t{dq}\n")
                else:
                    fout.write(f"{size}\n")
            fout.write("\n")

    print(f"Completed: generated {chain_id} chains")

if __name__ == '__main__':
    paf_to_chain(sys.argv[1], sys.argv[2])