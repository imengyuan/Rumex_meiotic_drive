import sys
import math

# calculating allele frequencies of different pools (including polarization)
# in sliding window (window size and step size 200kb)
# python3 freq_window.py scaffold.vcf scaffold_name start_pos
# version 1.0, input is from filt.vcf.py, after taking out sites of specific scaffold

f_in = open(sys.argv[1], "r")
out_file = sys.argv[1] + ".nobias.out"
f_out = open(out_file, "w")
# header line: scaffold name, window start, window end, avg freq in each window of 9 pools
print('scaffold\tstart_pos\tend_pos\tpool1', file=f_out)

# assign window size and step size
window = 200000  # 200kb
step = 200000
scaffold = sys.argv[2]
start_pos = int(sys.argv[3])
# initialize some variables
n = 1  # window NO.
result_list = []
# temporary variables
window_ref = [0]*9
window_alt = [0]*9
uncle_ref = [0]*9
uncle_alt = [0]*9
window_freq = [0]*9
line = f_in.readline()

while line:
    record = line.lstrip().split()
    if int(record[1]) < (n-1)*step + start_pos + window - 1:  # record[1] is the pos
        # calculate total read count in ref/alt for each pool
        for i in range(1, 2):  # for each pool(total9->6)
            window_ref[i] += int(record[2*i])
            window_alt[i] += int(record[2*i+1])
            uncle_ref[i] += int(record[5])
            uncle_alt[i] += int(record[6])
    else:
        start = (n-1)*step + start_pos
        end = (n-1)*step + start_pos + window - 1  # end = start + window - 1
        to_add = [scaffold, str(start), str(end)]
        # calculate allele freq
        for i in range(1, 2):

            # calculating mapping bias, constant across pools
            freq_R = uncle_ref[i] / float(uncle_ref[i] + uncle_alt[i])
            freq_N = 1 - freq_R
            bias = math.sqrt(freq_R * (1 - freq_N) / float(freq_N * (1 - freq_R)))

            # calculating freq of this window
            # notice the genome assembly is TX, we want Freq of NC alleles, so we calculate FreqN instead of FreR
            if window_ref[i] * bias + window_alt[i] != 0:
                # window_freq[i] = window_ref[i] * bias / float(window_ref[i] * bias + window_alt[i]) # NC is ref
                window_freq[i] = window_alt[i] / float(window_ref[i] * bias + window_alt[i])

            else:
                window_freq[i] = 0
                print(record[1], i)
            to_add.append(str(window_freq[i]))
        result_list.append(to_add)
        n += 1
        window_ref = [0] * 9  # initialize temporary variables
        window_alt = [0] * 9
        uncle_ref = [0] * 9
        uncle_alt = [0] * 9
        for i in range(1, 2):  # for each pool(total9->6)
            window_ref[i] += int(record[2*i])
            window_alt[i] += int(record[2*i+1])
            uncle_ref[i] += int(record[5])
            uncle_alt[i] += int(record[6])
    line = f_in.readline()

start = (n-1)*step + start_pos
end = start + window - 1
to_add = [scaffold, str(start), str(end)]
# calculate allele freq
for i in range(1, 2):
    # calculating mapping bias
    freq_R = uncle_ref[i] / float(uncle_ref[i] + uncle_alt[i])
    freq_N = 1 - freq_R
    bias = math.sqrt(freq_R * (1 - freq_N) / float(freq_N * (1 - freq_R)))
    # calculating freq of this window
    if window_ref[i] * bias + window_alt[i] != 0:
        # window_freq[i] = window_ref[i] * bias / float(window_ref[i] * bias + window_alt[i]) # NC is ref
        window_freq[i] = window_alt[i] / float(window_ref[i] * bias + window_alt[i])

    else:
        window_freq[i] = 0
        print(record[1], i)
    to_add.append(str(window_freq[i]))
result_list.append(to_add)

for i in result_list:
    print("\t".join(i), file=f_out)
f_in.close()
f_out.close()





