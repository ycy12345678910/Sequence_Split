import re

def find_main_modify(modify_list):
    depth,main_modify = 0,''
    if modify_list[0]:
        for i in range(len(modify_list[0]) - 1, -1, -1):  # 逆序遍历
            char = modify_list[0][i]
            if char in ['(','[','{']:
                    depth += 1
            elif char in [')',']','}']:
                depth -= 1
            if depth == 0:
                main_modify = modify_list[0][i:].strip('()').strip('[]').split(',') #去掉开头结尾的()[],按，拆分
                modify_list[0] = modify_list[0][:i]
                break
    return main_modify,modify_list

def sequence_mutation(main_modify,sequence_fundamental,AA_start_index):
    insert,AA_Symbol = '',''
    for chr in main_modify:
        #仅为数字时，为前一个氨基酸symbol
        if chr.isdigit():
            sequence_fundamental[int(chr)- AA_start_index] = AA_Symbol
        elif chr:
            #字符最末尾的数字为AA编号，其余为氨基酸symbol
            AA_Symbol = re.sub(r'\d+$', '', chr)
            match = re.search(r'(\d+)$', chr)[0]
            if 'Ace' in AA_Symbol:
                insert = 'Ace'
                AA_Symbol = AA_Symbol.replace('Ace-','')
            sequence_fundamental[int(match)- AA_start_index] = AA_Symbol
    return sequence_fundamental,insert

def sequence_split_modify(input_str):
    matches = re.search(r"[0-9a-zA-Z\-]+\((\d)+\-(\d+)\)",input_str)
    AA_start_index = int(matches.group(1))
    # AA_end_index = int(matches.group(2))
    fundamental_name = matches.group()
    modify_list = input_str.split(fundamental_name)  #按fundamental_name拆分字符串
    main_modify,modify_list = find_main_modify(modify_list)   #查找突变氨基酸的字符串
    #查找fundamental_name的序列，剪切，这里直接给定
    # sequence_fundamental = input()
    sequence_fundamental = ['Gly', 'Ser', 'Ser(n-octanoyl)', 'Phe', 'Leu', 'Ser', 'Pro', 'Glu', 'His', 'Gln', 'Arg', 'Val', 'Gln', 'Gln', 'Arg', 'Lys', 'Glu', 'Ser', 'Lys', 'Lys', 'Pro', 'Pro', 'Ala', 'Lys', 'Leu', 'Gln', 'Pro', 'Arg']
    sequence_modify,insert = sequence_mutation(main_modify,sequence_fundamental,AA_start_index)
    if modify_list[0]:
        #查找第一个数字，为修饰氨基酸
        modify_num = int(re.search('\d+',modify_list[0]).group())- AA_start_index
        sequence_modify[modify_num] += '-'+re.sub(r"-$", "", modify_list[0])
    if modify_list[1]:
        sequence_modify.append(re.sub(r"^-", "", modify_list[1]))
    
    #如果有封端，插入封端
    if insert:
        sequence_modify = [insert] + sequence_modify
        insert= ''
    # 过滤掉包含"des"的元素
    sequence_modify = [item for item in sequence_modify if "des-" not in item]
    return  sequence_modify