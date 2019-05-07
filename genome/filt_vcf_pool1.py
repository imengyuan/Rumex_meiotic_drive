import sys
import re

# input: output.vcf

f_in = open(sys.argv[1], "r")
out_file = sys.argv[1] + ".flit_vcf"
f_out = open(out_file, "w")
cnt_uncle = 0
cnt_DP_filt = 0
cnt_all = 0

content = f_in.readlines()
for line in content:  # loop through each line of the vcf
    if re.match('^#(.*)', line):  # skip headers(lines start with #)
        continue
    result_list = []  # a list keeping outputs of each line
    vcf_line = line.rstrip().split()  # eg. ['chrom_1','pos_1','.','G','A',...]
    uncle = vcf_line[16]
    S3 = vcf_line[18]
    S4 = vcf_line[19]
    S5 = vcf_line[20]
    S6 = vcf_line[21]

    # pool_1 = vcf_line[9]
    # new cut off of DP for pool 1
    pool_1_match = re.match('^.\/.\:(.+)\:(.+)\,(.+)$', vcf_line[9])
    pool_1_DP = int(pool_1_match.group(1))
    if (pool_1_DP >= 8) and (pool_1_DP <= 17):
        cnt_DP_filt += 1

        # filter for two genotype patterns
        if (uncle[0:3] == '0/1' and S3[0:3] == '0/0' and S4[0:3] == S5[0:3] == S6[0:3] == '1/1') or \
                (uncle[0:3] == '0/1' and S3[0:3] == '1/1' and S4[0:3] == S5[0:3] == S6[0:3] == '0/0'):
            # skip site where uncle's depth for heterozygous alleles is low
            # can be deleted
            uncle_match = re.match('^.\/.\:(.+)\:(.+)\,(.+)$', uncle)
            uncle_cnt_ref = int(uncle_match.group(2))
            uncle_cnt_alt = int(uncle_match.group(3))
            if uncle_cnt_alt == 0 or uncle_cnt_ref == 0:
                cnt_uncle += 1
                continue
            if uncle[0:3] == '0/1' and S3[0:3] == '0/0' and S4[0:3] == S5[0:3] == S6[0:3] == '1/1':
                uncle_cnt_ref, uncle_cnt_alt = uncle_cnt_alt, uncle_cnt_ref
            # first two column
            result_list.append(vcf_line[0])  # CHROM
            result_list.append(vcf_line[1])  # POS
            # output the genotype data, and polarize
            i = 9  # pool1
            if re.match('0\/0', vcf_line[i]) \
                    or re.match('0\/1', vcf_line[i]) or re.match('1\/1', vcf_line[i]):
                # extract two read counts in AD: cnt_ref, cnt_alt
                # GT:DP:AD        0/0:10:10,.
                AD_match = re.match('^.\/.\:(.+)\:(.+)\,(.+)$', vcf_line[i])
                cnt_DP = AD_match.group(1)
                cnt_ref = AD_match.group(2)
                cnt_alt = AD_match.group(3)
                if cnt_alt == ".":
                    cnt_alt = str(0)
                # Polarization
                if uncle[0:3] == '0/1' and S3[0:3] == '0/0' and S4[0:3] == S5[0:3] == S6[0:3] == '1/1':
                    cnt_ref, cnt_alt = cnt_alt, cnt_ref
                result_list.append(cnt_ref)
                result_list.append(cnt_alt)
                result_list.append(cnt_DP)
            else:  # cases of './.' allele freq = 0
                result_list.append("miss")  # no miss in pileup
                result_list.append("miss")

            result_list.append(str(uncle_cnt_ref))
            result_list.append(str(uncle_cnt_alt))
            cnt_all += 1
            print("\t".join(result_list), file=f_out)

f_in.close()
f_out.close()
print("# of filtered sites with uncle not true heterozygous: ", cnt_uncle)
print("# of filtered sites of DP>=8 and <=17: ", cnt_DP_filt)
print("# of all sites counted: ", cnt_all)

