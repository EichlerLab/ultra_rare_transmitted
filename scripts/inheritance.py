#!/usr/bin/env python
import argparse
import gzip
import sys

# USAGE STATEMENT AND HELP
parser = argparse.ArgumentParser(description='Determine parent of origin for informative sites in quad families.  Requires a VCF file and PED igree file as input.  Assumes samples have a sample ID structure of ID.fa, ID.mo, ID.p1, ID.s1.')
parser.add_argument('infile', help='input - tab-delimited VCF file') # POSITIONAL ARGUMENTS
parser.add_argument('ped', help = 'input - pedigree file') # POSITIONAL ARGUMENT
parser.add_argument('famid', help = 'input - family id') # POSITIONAL ARGUMENT
args = parser.parse_args() # PARSE ARGUMENTS

# INITALIZE INHERITANCE DICTIONARY 
geno_dict = {}

for f1 in ['1', '0']:
    for f2 in ['1', '0']:
        for m1 in ['1', '0']:
            for m2 in ['1', '0']:
                fa = '%(a1)s/%(a2)s' % {"a1": f1, "a2": f2}
                mo = '%(a1)s/%(a2)s' % {"a1": m1, "a2": m2}
                child = list(set(['%(a1)s/%(a2)s' % {"a1": x, "a2": y} for x,y in [(x, y) for x in [f1, f2] for y in [m1, m2]]] + ['%(a1)s/%(a2)s' % {"a1": x, "a2": y} for x,y in [(x, y) for x in [m1, m2] for y in [f1, f2]]]))
                geno_dict[fa + ':' + mo] = child

### FUNCTIONS ###
def check_mendelian(fa, mo, child, genos, sex, chrom):
    child = child.split(':')[0]
    if fa in ['0/0', '0/1', '1/0', '1/1'] and mo in ['0/0', '0/1', '1/0', '1/1'] and child in ['0/0', '0/1', '1/0', '1/1']:
        sex = sex.lower()[0]
        # ESTABLISH LIST OF POSSIBLE MENDELIAN GENOTYPES ACCORDING TO PARENTS AND OFFSPRING SEX
        if not any([missing_gt(fa), missing_gt(mo), missing_gt(child)]) and chrom != 'Y': # AUTOSOMES AND X CHROMOSOME
            ok = genos[fa + ':' + mo] # AUTOSOMES  
            if chrom in ['X', 'chrX']:
                fa = male_sex_chr(fa)
                if sex == 'm': # MALE X CHROMOSOME
                    child = male_sex_chr(child)
                    ok = mo.split('/')
                elif sex == 'f': # FEMALE X CHROMOSOME
                    ok = [a1 + '/' + a2 for a1 in fa.split('/') for a2 in mo.split('/') if a1 in ['1', '0']]
                    if '1/0' in ok:
                        ok.append('0/1')
                    elif '0/1' in ok:
                        ok.append('1/0')
        elif not any([missing_gt(fa), missing_gt(child)]) and chrom not in ['Y', 'chrY']: # Y CHROMOSOME   
            fa = male_sex_chr(fa)
            mo = './.'
            if sex == 'm': # MALE Y CHROMOSOME
                child = male_sex_chr(child)
                ok = fa[0]
            elif sex == 'f': # FEMALE Y CHROMOSOM
                child = '0/0'
                ok = '0/0'
        else:
            return(False)
        # CHECK CHILD GENOTYPES AGAINST LIST OF MENDELIAN GENOTYPES
        if child in ok: 
            return(True)
        else:
            return(False)
    else:
        return(False)

def male_sex_chr(gt):
    alleles = gt.split('/')
    if '1' in alleles:
        return('1')
    else:
        return('0')

def missing_gt(gt):
    if '.' in gt.split('/') or gt == '.':
        return(True)
    else:
        return(False)

# parental genotypes are formatted as father:mother
inh_dict = {'m':{'autosome':{'0/1:0/0': 'fa', '1/0:0/0': 'fa', '0/0:0/1': 'mo', '0/0:1/0': 'mo'}, 'X':{'0:0/1': 'mo', '0:1/0': 'mo', '1:1/0': 'mo'}, 'Y':{'1:0/0': 'fa'}},
'f':{'autosome':{'0/1:0/0': 'fa', '1/0:0/0': 'fa', '0/0:0/1': 'mo', '0/0:1/0': 'mo'}, 'X':{'1:0/0': 'fa', '0:1/0': 'mo', '0:0/1': 'mo'}}} 

def inheritance(fa, mo, child, sex, chrom):
    child = child.split(':')[0]
    if fa in ['0/0', '0/1', '1/0', '1/1'] and mo in ['0/0', '0/1', '1/0', '1/1'] and child in ['0/0', '0/1', '1/0', '1/1']:
        sex = sex.lower()[0]
        inh = 'not_inform'
        if chrom in ['X', 'chrX']:
            fa = male_sex_chr(fa)
            chrom = 'X'
            if sex == 'm':
                child = male_sex_chr(child)
        elif chrom in ['Y', 'chrY']:
            fa = male_sex_chr(fa)
            mo = '0/0'
            chrom = 'Y'
            if sex == 'm':
                child = male_sex_chr(child)
            elif sex == 'f':
                child = '0/0'
        else:
            chrom = 'autosome'        
        parents = fa + ':' + mo
        if '1' in child and parents in inh_dict[sex][chrom].keys():
            inh = inh_dict[sex][chrom][parents]
        else:
            inh = 'ref'
        return(inh)
    else:
        return('not_inform')

# POPULATE PHENOTYPE (GENDER/SEX) DICTIONARY
famid = args.famid
ped = {}
with open(args.ped) as f:
    for line in f:
        if line.split(' ')[0] == famid:
            dat = line.split(' ')
            ped[dat[1]] = dat[5].rstrip()

# VCF ANALYSIS
with gzip.open(args.infile) as f:
    for line in f:
        if line[0:2] == '##':
            sys.stdout.write(line)
        elif line[0] == '#':
            header = line.rstrip().split('\t')
            sys.stdout.write('##source=inheritance.py ' + args.infile + ' ' + args.ped + ' ' + args.famid + '\n')
            sys.stdout.write('##INFO=<ID=INHER,Number=1,Type=Character,Description="Inheritance pattern">' + '\n')
            sys.stdout.write('##INFO=<ID=MENDEL,Number=1,Type=Character,Description="Mendelian inheritance check">' + '\n')
            sys.stdout.write(line)
            
            #EXCEPTION HANDLING FOR FAMILIES MISSING PARENTS
            if 'fa' not in '\t'.join(header):
                sys.exit("No father present. Inheritance cannot be mapped.")
            if 'mo'  not in '\t'.join(header):
                sys.exit("No mother present. Inheritance cannot be mapped.")
            
            fa_idx = int([header.index(c) for c in header if famid+'.fa' == c][0])
            mo_idx = int([header.index(c) for c in header if famid+'.mo' in c][0])
            p_idx = map(int, [header.index(c) for c in header if famid+'.p' in c])
            s_idx = map(int, [header.index(c) for c in header if famid+'.s' in c])
            info = header.index('INFO')
        else:
            dat = line.rstrip().split('\t')
            fa = dat[fa_idx].split(':')[0]
            mo = dat[mo_idx].split(':')[0]
            mendel_p = all([check_mendelian(fa, mo, dat[p], geno_dict, ped[header[p]], dat[0]) for p in p_idx])
            inh_p = [inheritance(fa, mo, dat[p], ped[header[p]], dat[0]) for p in p_idx]
            if len(s_idx) > 0 and all(header[s] in ped.keys() for s in s_idx) and not any([missing_gt(dat[s]) for s in s_idx]):
                mendel_s = all([check_mendelian(fa, mo, dat[s], geno_dict, ped[header[s]], dat[0]) for s in s_idx])
                inh_s = [inheritance(fa, mo, dat[s], ped[header[s]], dat[0]) for s in s_idx]
            else:
                mendel_s = True
                inh_s = ['NA']
            mendel = all([mendel_p, mendel_s])
            inh = 'not_inform'
            for i in ['fa', 'mo']: # adapt for multiplex
                inh_list = []
                for p in range(0, len(inh_p)):
                    if inh_p[p] == i:
                        inh_list.append(header[p_idx[p]].split('.')[1])
                if len(s_idx) > 0:
                    for s in range(0, len(inh_s)):
                        if inh_s[s] == i:
                            inh_list.append(header[s_idx[s]].split('.')[1])
                if len(inh_list) > 0:
                    inh = i + '-' + ','.join(inh_list)
            if mendel and inh != 'not_inform':
                dat[info] = dat[info] + ';MENDEL=' + str(mendel) + ';INHER=' + inh
            else:
                dat[info] = dat[info] + ';MENDEL=' + str(mendel)
            # OUTPUT TO VCF
            sys.stdout.write('\t'.join(map(str, dat)) + '\n')       
