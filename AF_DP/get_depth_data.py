import sys
import re
# usage: python3 get_depth_data.py XX.vcf.table 1
# in this test, for pool sum2 the number should be 5
# bash loop as following
# for i in {1..5..1}; do python3 get_depth_data.py XX.vcf.table $i;done

f_in = open(sys.argv[1], "r")
pool_i = int(sys.argv[2])
out_file = sys.argv[1] + ".pool_" + sys.argv[2] + ".data"
f_out = open(out_file, "w")
# print header line of the output
print('DP\tN_SNP\tAvg_freq', file=f_out)
depths = []  #
freqs = []

content = f_in.readlines()
for line in content:
    if re.match('^#(.*)', line):  # skip header
        continue
    line = line.split()
    freqs.append(float(line[2*pool_i]))  # freq of pool i
    depths.append(int(line[2*pool_i+1]))  # depth of pool i


depth_values = set(depths)  # keep the unrepeated values of depths

depth_freq = [([0 for j in range(2)]) for i in range(len(depths))]  # so there're no index errors
for i in range(len(depths)):
    depth_freq[i][0] = depths[i]
    depth_freq[i][1] = freqs[i]
# print(depth_freq)

for DP in depth_values:
    N_SNP = depths.count(DP)  # the number of SNPs of this depth
    # calculating total_freq
    total_freq = 0
    for i in range(len(depths)):  # for all site in input table
        if depth_freq[i][0] == DP:
            total_freq += depth_freq[i][1]
    avg_freq = total_freq / N_SNP
    print(DP, '\t', N_SNP, '\t', avg_freq, file=f_out)

f_in.close()
f_out.close()

