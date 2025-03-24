import re
from SequenceSplit_Triple import sequence_split
amino_acid_map = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro","S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"}

def single2triple(single_S,sequence):
    for string in single_S:
        if string.isalpha():
            temp = amino_acid_map[string.upper()]
        elif re.search(r'\d|yl', string):  #有数字或者yl，则只转换第一个字符
            temp = amino_acid_map[string[0]] + string[1:]
        else:
            temp = "".join([amino_acid_map[char] if char.isalpha() else char for char in string])
        sequence.append(temp)
    return sequence

def sequence_split_Single(input_str):
    sequence = []
    split_str = sequence_split(input_str)  #按'-'切分
    for string in split_str:
        string_outside_parentheses = re.sub(r'[\(\[\{].*?[\)\]\}]', '', string) #判断括号外的字符串同时包含大写和小写/数字,或者都在（）里，则为修饰        
        if re.search(r'[A-Z]', string_outside_parentheses) and re.search(r'[a-z0-9]', string_outside_parentheses) or not string_outside_parentheses:
            sequence.append(string)
        else:
            # 使用正则表达式匹配单个字符或包含括号的部分  Y-Aib-EGT-aMeF(2F)-TSDYSI-aMeL-LDEK((2-[2-(2-Amino-ethoxy)-ethoxy]-acetyl)2-(yGlu)-CO-(CH2)18-COzH)AQ-Aib-EFI-(D-Glu)-YLIEGGPSSGAPPPS-NH2
            sequence_cut,AA_Symbol = [],''
            depth = 0  # 用于追踪括号嵌套层数
            for index,char in enumerate(string):
                AA_Symbol += char
                if char in ['(','[','{']:
                    depth += 1
                elif char in [')',']','}']:
                    depth -= 1
                if depth == 0 and (index +1 < len(string) and string[index+1] not in ['(','{','[' ]) or index+1  == len(string) :
                    sequence_cut.append(AA_Symbol)
                    AA_Symbol = ''
            sequence = single2triple(sequence_cut,sequence)
    return sequence
