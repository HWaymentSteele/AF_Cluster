from glob import glob

#RfaH

possible_msas = glob('data_sep2022/04_OtherFoldswitchers/00_RfaH/msas/*.a3m')

lst=[]
for m in possible_msas:
    lins = open(m,'r').readlines()
    n_seqs = len(lins)/2 - 1
    if n_seqs>=10 and 'U' not in m:
        lst.append(m.split('_')[-1].replace('.a3m',''))
print(len(lst))

existing=['049','056','021','014','045']+['001','002','005','024','025']

print(sorted([x for x in lst if x not in existing]))

#Mad2

possible_msas = glob('data_sep2022/04_OtherFoldswitchers/01_Mad2/1S2H_msas/*.a3m')

lst=[]
for m in possible_msas:
    lins = open(m,'r').readlines()
    n_seqs = len(lins)/2 - 1
    if n_seqs>=10 and 'U' not in m:
        lst.append(m.split('_')[-1].replace('.a3m',''))
print(len(lst))

existing=['000','046','008']+['020','055','083']

print(sorted([x for x in lst if x not in existing]))
print(len(sorted([x for x in lst if x not in existing])))
