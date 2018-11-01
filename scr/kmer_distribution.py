#!/Users/ashleyli/anaconda3/envs/py3/bin/python

import sys
#deal with interleaved seq and get seq from a genome fasta file
def get_seq(input_file):
    seq_list = []
    with open(input_file , 'r') as fin:
        seq = ''
        for line in fin:
            line = line.rstrip()
            if line.startswith('>'):
                if seq != '':
                    seq_list.append(seq)
                    seq = ''
            else:
                seq += line
        seq_list.append(seq)
    return seq_list
   
#get kmer of each sequence
def get_kmer(seq, kmer_size):
    if kmer_size >= len(seq):
        print("\nError:kmer size greater than sequence length")
        sys.exit()
    i = 0
    kmer_list = []
    while i + kmer_size <= len(seq):
        kmer_list.append(seq[i:i+kmer_size])
        i += 1
    return kmer_list

#calculate the abundance of each kmer in the kmer list   
def kmer_dis(kmer_list):
    kmer_dict = {}
    kmer_set = set(kmer_list)
    for kmer in kmer_set:
        kmer_num = kmer_list.count(kmer)
        kmer_dict[kmer] = kmer_num
    return kmer_dict

#calculate frequency of each kmer    
def cal_frq(kmer_dict):
    total = 0
    for num in kmer_dict.values():
        total += num
    for kmer,num in kmer_dict.items():
        frq = num/total
        kmer_dict[kmer] = frq
    return kmer_dict
        
        
        
        


    
    
    
        
    

    
    

                
                
