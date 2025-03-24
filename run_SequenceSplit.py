import pandas as pd
import re
from SequenceSplit_Triple import sequence_split_triple
from SequenceSplit_Single import sequence_split_Single
from SequenceSplit_Mutation import sequence_split_modify

def contains_amino_acid(s):
    amino_acids = r"Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val"
    s = re.sub(r'\(.+?\)','',s)
    match = re.search(amino_acids, s)
    return bool(match) 

def sequence_split(string):
    if string =='':
        return ''
    if set(string).issubset({"A","G","C",'T','U'}) :
        return "DNA/RNA"
    string = re.sub(r'^NH2-|-OH$|-COOH$', '',string).replace('Ac','Ace').replace('L-','').replace(' ','').replace("CONH2","NH2")
    string = re.sub(r'(?<!-)(NH2)$', r'-NH2', string) #非-NH2结尾的NH2变成-NH2
    string = re.sub(r'(?<=[a-z])-(?=[a-z])', '', string)  #去除前后都为小写字母的-，为换行连接
    #判断序列类型是单字母为框架/三字母为框架/序列加上突变为框架
    if re.search(r"\(\d+\-\d+\)", string):  #匹配（数字-数字），为框架多肽+突变
        sequence = sequence_split_modify(string)
        # sequence = []
        pass
    elif contains_amino_acid(string): #三字母
        #替换H-开头的
        string = re.sub(r'^H-', '',string)
        sequence= sequence_split_triple(string)
    else:   #单字母
        sequence= sequence_split_Single(string)
    return sequence

# file_name = '../WO2007041278-1'
# df = pd.read_csv(f'{file_name}.csv')
# df['symbol-split'] = df['COMPOUND'].apply(sequence_split)
# df = df[['COMPOUND','symbol-split']]
# df.to_csv(f'{file_name}-smiles.csv', index=False)
