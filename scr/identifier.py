#!/Users/ashleyli/anaconda3/envs/py3/bin/python

import sys
import argparse
from kmer_distribution import get_seq, get_kmer, kmer_dis, cal_frq
from scipy.stats import chisquare,ks_2samp

file_num = len(sys.argv) - 8

parser = argparse.ArgumentParser(description="sample usage: identifier.py -k 2 -g Pk.fasta Pb.fasta Pv.fasta Py.fasta -s test -t 2")
parser.add_argument("-k", "--kmer", help = "Kmer size", type = int)
parser.add_argument("-g", "--genome", help = "Reference fasta file" , nargs=file_num)
parser.add_argument("-s", "--sequence", help = "Undetermined sequence")
parser.add_argument("-t", "--test", help = "Statistical test [1:x2 test; 2:k-s test]",type = int,choices = [1, 2])
args = parser.parse_args()

file_name_list= args.genome
kmer_size = args.kmer

#get a dictionary which contains the kmer frequency distribution of each reference genome
#return kmer_fullref_dict = {genome1:{kmer1:0.12,kmer2:0.033,...},genome2:{...},...}
kmer_fullref_dict = {}
for i in range(file_num):
    kmer_ref_list = []
    input_ref_file = args.genome[i]
    for seq in get_seq(input_ref_file):
        kmer_ref_list.extend(get_kmer(seq, kmer_size))
    kmer_ref_dict = kmer_dis(kmer_ref_list)
    ref_frq_dict = cal_frq(kmer_ref_dict)
    kmer_fullref_dict[args.genome[i].replace('.fasta','')] = ref_frq_dict

#get a dictionary which contains the kmer frequency distribution of undetermined sequence
#return test_frq_dict = {kmer1:0.123,kmer2:...}
kmer_test_list = []
input_test_file = args.sequence
for seq in get_seq(input_test_file):
    kmer_test_list.extend(get_kmer(seq, kmer_size))
kmer_test_dict = kmer_dis(kmer_test_list)
test_frq_dict = cal_frq(kmer_test_dict)

#find same kmers in ref genome and sample sequence
#return same_kmer_dict = {genome1:[kmer1,kmer2,...], genome2:[....],...}
same_kmer_dict = {}
kmer_test = test_frq_dict.keys()
for each_ref in kmer_fullref_dict:
    kmer_ref = kmer_fullref_dict[each_ref].keys()
    same_kmer =  set(kmer_test) & set(kmer_ref)
    same_kmer_dict[each_ref] = same_kmer

#x2 test
def x2_test(ref_dict,test_dict):
    ref_list = []
    test_list = []
    for kmer in ref_dict:
        ref_list.append(ref_dict[kmer])
        test_list.append(test_dict[kmer])
    #to create two list which in the same kmer order
    result = chisquare(test_list, f_exp = ref_list)
    return result
    
#k-s test
def k_s(ref_dict,test_dict):
    ref_list = []
    test_list = []
    for kmer in ref_dict:
        ref_list.append(ref_dict[kmer])
        test_list.append(test_dict[kmer])
    result = ks_2samp(test_list, ref_list)
    return result

#do the statistics tests
#return result_dict = {genome1:result(statistic_value, p_value), genome2:...}
result_dict = {}
for ref_name,kmer_list in same_kmer_dict.items():
    ref_dict = {}
    test_dict = {}
    if len(kmer_list) == 0:
        print('\nError:One of the reference gemomes and test sequence share no same kmer, please change your reference genome or choose a smaller kmer size')
        sys.exit()
    for kmer in kmer_list:
        ref_dict[kmer] = kmer_fullref_dict[ref_name][kmer]
        test_dict[kmer] = test_frq_dict[kmer]
    #to create two dictionaries(reference & test) which share the same kmers 
    if args.test == 1:
        result_dict[ref_name] = x2_test(ref_dict, test_dict)
    if args.test == 2:
        result_dict[ref_name] = k_s(ref_dict, test_dict)

#print the result       
if  args.test == 1:
    print('\nChi-squared test result:')
    print('\nref_name\tx2-value\tP-value')
if args.test == 2:
    print('\nKolmogorov–Smirnov test result:')
    print('\nref_name\tD-value\tP-value')   
for ref_name in result_dict:
    print(ref_name + '\t' + str(result_dict[ref_name][0]) + '\t' + str(result_dict[ref_name][1]))

#select the reference genome based on the statistics value
result = min(result_dict.items(), key=lambda x: x[1][0])
print('\nThe sample sequence could belong to:')
print(result[0] + '\tstatistic=' + str(result[1][0]) + '\tp-value=' + str(result[1][1]))
