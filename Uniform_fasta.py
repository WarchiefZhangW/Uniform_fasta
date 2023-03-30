#!/usr/bin/python
import os, sys
import multiprocessing as mp

"""

#############################################################
#******************Uniform Fasta File Tool******************#
#***********************by ZHANG W**************************#
#******zhangww@big.ac.cn    Wechat&TEL:18210123493**********#
#***********Beijing Institute of Genomics, CAS, CN**********#
#***********************************************************#
#############################################################
"""

def read_fasta(file_name):#a generator read the fasta file
        name = ''
        seq = ''
        with open(file_name) as rawdata:
                
                for line in rawdata:
                        if line=='' or line=='\n':
                                continue
                        if line[0]=='>':
                                if name != '':
                                        yield name, seq
                                name = line.split()[0][1:]
                                seq = ''
                        if line[0]!='>':
                                seq += line.strip()
        yield name, seq

def uniform_data(seq):
    CB = {chr(i):i for i in range(32, 127)}
    NB = {'A':0, 'T':0, 'C':0, 'G':0, 'N':0, 'a':0, 't':0, 'c':0, 'g':0, 'n':0}
    read_id, read_name, read_seq = seq
    new_name = ''
    new_seq = ''
    for w in read_name:
        if w not in CB:
            continue
        new_name += w.upper()
    for w in read_seq:
        if w not in NB:
            continue
        new_seq += w.upper()
    return read_id, new_name, new_seq



def uniform_data_single():
    CB = {chr(i):i for i in range(32, 127)}
    NB = {'A':0, 'T':0, 'C':0, 'G':0, 'N':0, 'a':0, 't':0, 'c':0, 'g':0, 'n':0}
    uniformed_fasta = open('{}_uniform.fasta'.format(sys.argv[1].split('.')[0]),'w')
    for i, (read_name, read_seq) in enumerate(read_fasta(sys.argv[1])):
        if i % 100 ==0:
            print('\r', 'Uniformed {} reads...'.format(i), end = '')
        #print(i, read_name)
        new_seq = ''
        new_name = ''
        for w in read_name:
            if w not in CB:
                continue
            new_name += w.upper()
        for w in read_seq:
            if w not in NB:
                continue
            new_seq += w.upper()
        uniformed_fasta.writelines('>{0}\n{1}\n'.format(new_name, new_seq))


def uniform_data_mp(t):
    uniformed_fasta = open('{}_uniform.fasta'.format(sys.argv[1].split('.')[0]),'w')
    old_data = []
    for i,(name,seq) in enumerate(read_fasta(sys.argv[1])):
        #print(name, i)
        old_data.append((i, name, seq))
        if i % t*4 == 0 and i != 0:
            print('\r', 'Uniformed {} reads...'.format(i), end = '')
            p = mp.Pool(processes = t)
            new_data = p.map(uniform_data, old_data)
            p.close()
            p.join()
            new_data.sort(key = lambda x:x[0])
            new_data_write = ['>{0}\n{1}\n'.format(new_name, new_seq) for _, new_name, new_seq in new_data]
            uniformed_fasta.writelines(new_data_write)
            old_data = []
            new_data = []
    if old_data != []:
        p = mp.Pool(processes = t)
        new_data = p.map(uniform_data, old_data)
        new_data.sort(key = lambda x:x[0])
        p.close()
        p.join()
        new_data_write = ['>{0}\n{1}\n'.format(new_name, new_seq) for _, new_name, new_seq in new_data]
        uniformed_fasta.writelines(new_data_write)
    uniformed_fasta.close()

def main():
    usage = 'Usage:  ./Uniform_fasta.py {file_name} {threads}'
    if len(sys.argv) == 1:
        print(usage)
        sys.exit()
    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print(usage)
        sys.exit()
    try: t = int(sys.argv[2])
    except: t = 1
    print('Uniforming fasta sequences, using {} threads'.format(t))
    if t == 1:
        uniform_data_single()
    elif t > 1:
        uniform_data_mp(t)

if __name__ == '__main__':
    main()
