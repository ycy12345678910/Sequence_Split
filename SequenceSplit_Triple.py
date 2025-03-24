import re

def modify_cut(str_cut):
    str_cut_modify,AA_Symbol = [],''
    for char in str_cut:
        # if len(char) < 3 or not re.search(r'[A-Z]', char): #避免D-,4-,4F-benzoyl-Arg等修饰的剪切
        char_out_parentheses = re.sub(r'[\(\[\{].*[\)\]\}]', '', char) 
        #避免D-,4-,4F-benzoyl-Arg,"D,L-(a-Methyl)Trp"等修饰的剪切
        if len(char) < 3 or 'yl' in char_out_parentheses or not(re.search(r'[A-Z]', char_out_parentheses) and re.search(r'[a-z0-9]', char_out_parentheses)) : 
            AA_Symbol += f'{char}-'
        else:
            AA_Symbol += char
            str_cut_modify.append(AA_Symbol)
            AA_Symbol = ''
    return str_cut_modify

def sequence_split(input_str):
    sequence_cut,AA_Symbol = [],''
    #对字符串进行切分，仅匹配括号外的 "-"
    depth = 0  # 用于追踪括号嵌套层数
    for char in input_str:
        if char in ['(','[','{']:
            depth += 1
        elif char in [')',']','}']:
            depth -= 1
        if char == '-' and depth == 0 :  # 仅在括号外部才切分
            sequence_cut.append(AA_Symbol)  
            AA_Symbol = ''  
        else:
            AA_Symbol += char
    sequence_cut.append(AA_Symbol) #添加最后一个氨基酸
    return sequence_cut

def sequence_split_triple(input_str):
    sequence_cut = sequence_split(input_str)
    sequence_cut = modify_cut(sequence_cut)
    AA_Symbol,sequence = '',[]
    #对切分的字符串整理为序列，以 ”（ 开头，包含字母和数字，且以 ）数字 结尾，这种情况数字表示重复次数
    for AA_Symbol in sequence_cut:
        match = re.search(r'^\(([a-zA-Z0-9]+)\)(\d+)$', AA_Symbol)
        if match:
            repeat_num = int(match.group(2))
            AA_Symbol = match.group(1)
            sequence.extend([AA_Symbol] * repeat_num)
        else:
            sequence.append(AA_Symbol)
    return sequence
