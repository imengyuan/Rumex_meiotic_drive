import sys
import re

# A script for filtering a combined vcf based on some rule of genotypes
# and calculate allele frequencies for each sample
# Usage: python3 get_freq.py test.combined.vcf
# Input file: test.combined.vcf
# Output file: test.combined.vcf.table
# pool sum2-> combined bam, pool sum-> individual bam's counts sum

f_in = open(sys.argv[1], "r")
out_file = sys.argv[1] + ".table"
f_out = open(out_file, "w")
# print header line of the output
print('#CHROM\tPOS\tpool_1\tpool_2\tpool_3\tpool_4\tpool_5\tpool_6\tpool_7\tpool_combine\tpool_sum', file=f_out)
cnt_stars = 0  # number of star sites
cnt_uncle = 0
cnt_all = 0
cnt_DP = 0

content = f_in.readlines()
for line in content:  # loop through each line of the vcf
    if re.match('^#(.*)', line):  # skip headers(lines start with #)
        continue
    result_list = []  # a list keeping outputs of each line
    vcf_line = line.split()  # eg. ['chrom_1','pos_1','.','G','A',...]
    uncle = vcf_line[16]
    S3 = vcf_line[18]
    S4 = vcf_line[17]
    S5 = vcf_line[19]
    S6 = vcf_line[20]
    # skip sites with * in alt (likely indels)
    if vcf_line[4] == "*":  # alt
        cnt_stars += 1  # no indels for mpileup, can be deleted
        continue
    # identify ancestry informative sites (two patterns included)
    if (uncle[0:3] == '0/1' and S3[0:3] == '0/0' and S4[0:3] == S5[0:3] == S6[0:3] == '1/1') or \
            (uncle[0:3] == '0/1' and S3[0:3] == '1/1' and S4[0:3] == S5[0:3] == S6[0:3] == '0/0'):
        # filter sites with DP>2*mean DP
        DP_pools = 0
        for i in range(9, 16):  # for 7 pools no pool combine
            match = re.match('^.\/.\:(.*)\,(.*?)\:.*', vcf_line[i])
            DP_pools += int(match.group(1)) + int(match.group(2))
        mean_DP_pools = DP_pools / 8.0000
        if mean_DP_pools > 245.9403:
            cnt_DP += 1
            continue

        sum_ref = 0
        sum_alt = 0
        uncle_AD_match = re.match('^.\/.\:(.*)\,(.*?)\:.*', uncle)
        uncle_cnt_ref = int(uncle_AD_match.group(1))
        uncle_cnt_alt = int(uncle_AD_match.group(2))
        if uncle_cnt_alt == 0 or uncle_cnt_ref == 0:
            cnt_uncle += 1
            continue
        result_list.append(vcf_line[0])  # CHROM
        result_list.append(vcf_line[1])  # POS
        # calculating allele frequencies of different pools (including polarization)
        for i in range(9, 22):  # for 7 samples and pool sum2
            if i in range(16, 21):  # skip the rest
                continue
            if re.match('0\/0', vcf_line[i]) \
                    or re.match('0\/1', vcf_line[i]) or re.match('1\/1', vcf_line[i]):
                # extract two read counts in AD: cnt_ref, cnt_alt
                AD_match = re.match('^.\/.\:(.*)\,(.*?)\:.*', vcf_line[i])
                cnt_ref = int(AD_match.group(1))
                cnt_alt = int(AD_match.group(2))
                # Polarization: inverse frequency values 0.75->0.25
                if uncle[0:3] == '0/1' and S3[0:3] == '1/1' and S4[0:3] == S5[0:3] == S6[0:3] == '0/0':
                    cnt_ref, cnt_alt = cnt_alt, cnt_ref
                # calculate allele frequency
                if cnt_ref + cnt_alt != 0:
                    allele_freq = cnt_ref / float(cnt_ref + cnt_alt)
                else:
                    print("0 ", i, " ", line.rstrip())
                    allele_freq = 0  # need a check
                result_list.append(str(allele_freq))
                # do not add counts from pool sum2!
                if i == 22:
                    cnt_ref, cnt_alt = 0, 0
                sum_ref += cnt_ref
                sum_alt += cnt_alt
            else:  # cases of './.' allele freq = 0 (need a check)
                print("./. ", i, " ", line.rstrip())
                result_list.append(str(0))
        if sum_ref + sum_alt != 0:
            sum_allele_freq = sum_ref / float(sum_ref + sum_alt)
            result_list.append(str(sum_allele_freq))
        else:
            print("sum0 ", line.rstrip())
            result_list.append(str(0))
        # if result_list is not empty, write the results to output
        cnt_all += 1
        print("\t".join(result_list), file=f_out)
f_in.close()
f_out.close()
print("# of filtered sites with * in alt: ", cnt_stars)
print("# of filtered sites with uncle not true heterozygous: ", cnt_uncle)
print("# of filtered sites with DP > 2*mean_DP: ", cnt_DP)
print("# of all sites counted: ", cnt_all)

