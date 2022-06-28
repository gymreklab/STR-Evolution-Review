# imports
import pandas as pd
import numpy as np
from datetime import datetime
import sys

if len(sys.argv) < 3:
    print("Usage: python filter_TRF_nopandas.py in.bed out.bed")
    sys.exit()
in_bed_file = sys.argv[1]
out_bed_file = sys.argv[2]
# in_bed_file = 'trf_out/repeats_chr1_2.txt'
# out_bed_file = 'trf_out/repeats_chr1_2_filt_local.txt'

def get_col_as_list(mylist, col):
    vals = []
    for i in mylist:
        vals.append(getattr(i, col))
    return(vals)

# group_by_col = 'locus'
def keep_smallest_motif(line_list, group_by_col):

    # get list of group ids for every row
    gids = get_col_as_list(line_list, group_by_col)
    idx = list(range(len(gids)))

    # find rows for which ids are duplicated
    gids, idx = zip(*sorted(zip(gids, idx)))
    is_dup = []
    for i in range(1,len(gids)-1):
        if gids[i] == gids[i+1] or gids[i] == gids[i-1]:
            is_dup.append(True)
        else:
            is_dup.append(False)
    is_dup = [gids[0] == gids[1]] + is_dup
    is_dup = is_dup + [gids[-2] == gids[-1]]
    idx, is_dup = zip(*sorted(zip(idx, is_dup)))

    # select data lines where group id is part of the set of duplicates
    dup_list = [line_list[i] for i,d in enumerate(is_dup) if d]

    # sort dup_list by the grouping variable in case it is not sorted
    dup_gids = get_col_as_list(dup_list, group_by_col)
    idx = list(range(len(dup_gids)))
    if len(dup_gids)==0:
        return (line_list)
    # rearrange dup_list according to the new index
    sort_dup_gids, sort_idx = zip(*sorted(zip(dup_gids, idx)))
    dup_list = [dup_list[i] for i in sort_idx]

    # loop over each row
    # keep track of the motif and the associated rowid
    # once a new locus is encountered save the rowid associated with shorted motif to keep
    # keep multiple in case of ties
    i = 0
    keep = []
    discard = []
    loc = getattr(dup_list[i], group_by_col)
    mot = [dup_list[i].motif]
    row = [dup_list[i].rowid]
    while i < (len(dup_list)-1):
        next_loc = getattr(dup_list[i+1], group_by_col)
        next_mot = dup_list[i+1].motif
        next_row = dup_list[i+1].rowid
        if next_loc == loc:
            mot.append(next_mot)
            row.append(next_row)
        else:
            # which_min = np.argmin([len(m) for m in mot])
            # keep.append(row.pop(which_min))
            # discard = discard + row
            ml = np.array([len(m) for m in mot])
            keep = keep + np.array(row)[ml == ml.min()].tolist()
            discard = discard + np.array(row)[ml != ml.min()].tolist()
            mot = [next_mot]
            row = [next_row]
        loc = next_loc
        i += 1
    # finish out the last comparison
    ml = np.array([len(m) for m in mot])
    keep = keep + np.array(row)[ml == ml.min()].tolist()
    discard = discard + np.array(row)[ml != ml.min()].tolist()

    # dicard the lines with rowid in discard from the main list
    line_dict = dict()
    for l in line_list:
        line_dict[l.rowid] = l
    to_keep_main = set(line_dict.keys()) - set(discard)
    deduped = []
    for id in to_keep_main :
        deduped.append(line_dict[id])
    len(deduped)

    return(deduped)

# Check if the first and last K bps match, trim half motifs at the end.
def count_motif_tandem(repeat_str, motif):
    m = len(motif)
    c = 0
    s = 0
    p = []
    while s < len(repeat_str):
        #print(repeat_str[s:s + m])
        if (repeat_str[s:s + m] == motif):
            c = c + 1
            p.append(s)
            s = s + m
        else:
            s = s + 1
    spacing = [x-y for x,y in zip(p[1:len(p)], p[0:len(p)-1])]

    rep_runs = [1]
    for i in range(len(p)-1):
        if spacing[i] == m:
            rep_runs[len(rep_runs)-1] += 1
        else:
            rep_runs.append(1)

    if max(rep_runs) > 1:
        c_tandem = max(rep_runs)
    else:
        c_tandem = 0

    return c_tandem

# Find coumpund motifs: ATAT = (AT)*2
def is_compound_tandem(motif):
    l = len(motif)
    threshold = 1
    for i in range (1, int(l / 2) + 1):
        sub = motif[0:i]
        # print count_motif(motif, sub) * i, l * threshold
        if count_motif_tandem(motif, sub) * i == l * threshold :
            return True
    return False

def minimal_trim(rep, motif):
    mm = motif * 2
    ll = len(motif) * 2
    start_match = False
    end_match = False
    max_trim_len = len(rep)-ll
    for start_offset in range(max_trim_len + 1):
        if (start_offset + ll > len(rep)):
            return -1,-1
        if rep[start_offset: start_offset + ll] == mm:
            start_match = True
            break
    if start_match == False:
        return -1, -1
    for end_offset in list(reversed(range(ll, len(rep)+1))):
        if end_offset - ll < 0:
            return -1, -1
        if rep[end_offset - ll: end_offset] == mm:
            end_match = True
            break
    if end_match == False:
        return -1, -1
    return start_offset, end_offset

def expand_string(s, times):
    rem = np.mod(times, 1)
    whl = int(times - rem)
    expanded = s*whl + s[0:int(rem*len(s))]
    return(expanded)

# time script
startTime = datetime.now()

class data_line:
    def __init__(self, vals):
        # assign vars
        self.chr = vals[0]
        self.start = int(vals[1])
        self.end = int(vals[2])
        self.motif_len = len(vals[4])
        self.motif = vals[4]
        self.rep_string = vals[5]
        self.rowid = vals[6]

        # create new vars
        self.locus = self.chr + '_' + str(self.start) + '_' + str(self.end)
        self.start_id = self.chr + '_' + str(self.start)
        self.end_id = self.chr + '_' + str(self.end)

    def __str__(self):
        return('\n'.join([(str(k) + ': ' + str(v)) for k,v in self.__dict__.items()]))

    def copy(self):
        new_data_line = type(self)([self.chr, self.start, self.end, self.motif_len, self.motif, self.rep_string, self.rowid])
        return(new_data_line)

print(str.format('Reading data from {0}', in_bed_file))
# dictionary entry for each chromosome
chrs = dict()
# read in the repeats from TRF
with open(in_bed_file, 'r') as in_bed:
    # set up the first line list from the first line of the file
    line_list = []; lnum = 1
    line = in_bed.readline().strip().split('\t')
    curr_chr = line[0]
    line_list.append(data_line(line + [lnum]))
    while True:
        # read new line
        line = in_bed.readline()
        if not line:
            break

        line = line.strip().split('\t')
        lnum += 1

        # start a new chromosome if necessary and add previous to dictionary and reset line_list
        if line[0] != curr_chr:
            chrs[curr_chr] = line_list
            curr_chr = line[0]
            line_list = []

        # append the line to line_list
        line_list.append(data_line(line + [lnum]))
    # add last list
    chrs[curr_chr] = line_list

# time script
#loadTime = datetime.now()
#print(str.format('Time to load input data\t{0}', loadTime - startTime))

# process chromosome by chromosome
first_out = True
for chr, line_list in chrs.items():
    print(str.format('------------------\nProcessing chromosome: {0}\n------------------', chr))

    print(str.format('Number of repeats: {0}', len(line_list)))

    # first reduce any sets of repeats that share a common start and end to one
    # same locus: same chr, start and end
    # keep the shortest (simplest) motif
    if len(line_list) >1:
        line_list = keep_smallest_motif(line_list, 'locus')
    print(str.format('Keep smallest motif from identical loci; Survived: {0} repeats', len(line_list)))

    # second reduce any sets of repeats that share a common start but not end
    if len(line_list) >1:
        line_list = keep_smallest_motif(line_list, 'start_id')
    print(str.format('Keep smallest motif from identical start pos; Survived: {0} repeats', len(line_list)))

    # second reduce any sets of repeats that share a common end but not start
    if len(line_list) >1:
        line_list = keep_smallest_motif(line_list, 'end_id')
    print(str.format('Keep smallest motif from identical end pos; Survived: {0} repeats', len(line_list)))

    # time script
    # discardTime = datetime.now()
    # print(str.format('Time to reduce duplicatd repeats\t{0}', discardTime - loadTime))

    # define threshold for minimum number of repeats
    thresholds = {1:10, 2:5, 3:4, 4:3, 5:3, 6:3}

    print('Trimming repeats, filtering homopolymers and compound repeats')
    r = 0
    proc_line_list = []
    for line in line_list:
        # if np.mod(r, 10000) == 0:
        #     print(r)
        r += 1

        # remove homopolymers
 #       if line.motif_len == 1:
 #           continue

        # check if motif is compound
        if is_compound_tandem(line.motif):
            continue

        # find the minimally trimmed sequence two ways
        st, en = minimal_trim(line.rep_string, line.motif)

        # adjust the repeat string
        if st == -1 or en == -1:
            continue

        line.rep_string = line.rep_string[st:en]
        line.end = line.start + en - 1
        line.start = line.start + st
        ref_copy = (line.end - line.start + 1) / line.motif_len

        # threshold
        if line.motif_len in thresholds:
            thresh = thresholds[line.motif_len]
        else:
            thresh = 3

        # remove those motifs which don't meet the threshold
        if (ref_copy < thresh):
            continue

        # remove non perfect motifs
        if (expand_string(line.motif, ref_copy) != line.rep_string):
            continue

        # if passed everything add to proc list
        proc_line_list.append(line)

    print(str.format('Survived: {0} repeats', len(proc_line_list)))
    print('Removing duplicate repeats')
    #check if survived repeats equal to 0
    if len(proc_line_list) == 0:
        print('No repeats found ... skipping to next chrom')
        continue
    # sort proc_line_list by chr_start before overlap detection
    start_idx = get_col_as_list(proc_line_list, 'start')
    idx = list(range(len(start_idx)))
    
    # rearrange dup_list according to the new index
    sort_start_idx, sort_idx = zip(*sorted(zip(start_idx, idx)))
    proc_line_list = [proc_line_list[i] for i in sort_idx]

    # find all instances of duplicated repeats
    dup_line_idx = []
    for r in range(len(proc_line_list)-1):
        line = proc_line_list[r]
        next_line = proc_line_list[r+1]
        if line.start == next_line.start and line.end == next_line.end:
            dup_line_idx.append(r)

    # remove duplicated repeats
    to_keep = set(range(len(proc_line_list))) - set(dup_line_idx)
    dedup_list = []
    for i in to_keep:
        dedup_list.append(proc_line_list[i])
    

    print(str.format('Survived: {0} repeats', len(dedup_list)))
    print('Removing overlapping repeats')
    #check if survived repeats equal to 0
    
    if len(dedup_list) == 0:
        print('No repeats found ... skipping to next chrom')
        continue
    # find instance with overlap
    overlap_line_idx = []
    for r in range(len(dedup_list)-1):
        line = dedup_list[r]
        next_line = dedup_list[r+1]
        if line.start == next_line.start and line.end == next_line.end:
            dup_line_idx.append(r)

        # check that there isn't overlap between the repeats
        if next_line.start < line.end:
            # check that the bases used are similar, i.e. are these a pair of ambiguous repeats
            # by taking symmetric different of the sets of motif letters, i.e. in A or B, but not both
            if set(next_line.motif) ^ set(line.motif) == set():
                overlap_line_idx.append(r)

    # grab the next line from any line in overlap line index as well
    overlap_line_idx = overlap_line_idx + [i + 1 for i in overlap_line_idx]
    overlap_line_idx = list(set(overlap_line_idx))

    # remove overlap lines
    to_keep = set(range(len(dedup_list))) - set(overlap_line_idx)
    unoverlap_list = []
    for i in to_keep:
        unoverlap_list.append(dedup_list[i])

    print(str.format('Survived: {0} repeats', len(unoverlap_list)))
    print(str.format('Writing output file to: {}', out_bed_file))

    #
    if first_out:
        fotype = 'w'
        first_out = False
    else:
        fotype = 'a'

    # write the final file
    with open(out_bed_file, fotype) as outfile:
        for line in unoverlap_list:
            out_line = str.join('\t', [line.chr, str(line.start), str(line.end), str(line.motif_len), line.motif,
                                       line.rep_string]) + '\n'
            outfile.write(out_line)

# time script
endTime = datetime.now()
print(str.format('Done! Processing time was: {0}', endTime - startTime))
