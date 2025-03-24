import json
import re
import copy
import os,time
from Symbol2Smiles import symbol_to_smiles 
from run_SequenceSplit import sequence_split
symbol_dict = {'Apc': 'NC1(C(=O)O)CCNCC1',
               'D-1-Nal': 'N[C@@H](CC1=CC=CC2=CC=CC=C12)C(=O)O',
               'D-Trp': 'O=C(O)[C@H](N)Cc1c[nH]c2c1cccc2',
               '2-Thi': 'N[C@@H](CC1=CC=CS1)C(=O)O',
               'Inp': 'C1CNCCC1C(=O)O',
               'D-2-Nal': 'N[C@@H](CC1=CC2=CC=CC=C2C=C1)C(=O)O',
               'D-Bal': 'N[C@H](CC1=CSC2=CC=CC=C12)C(=O)O',
               'Ala': 'C[C@@H](C(=O)O)N',
               'Cys': 'C([C@@H](C(=O)O)N)S',
               'Asp': 'C([C@@H](C(=O)O)N)C(=O)O',
               'Glu': 'C(CC(=O)O)[C@@H](C(=O)O)N',
               'Gly': 'C(C(=O)O)N',
               'His': 'O=C(O)[C@@H](N)Cc1nc[nH]c1',
               'Ile': 'CC[C@H](C)[C@@H](C(=O)O)N',
               'Leu': 'CC(C)C[C@@H](C(=O)O)N',
               'Lys': 'NCCCC[C@H](N)C(=O)O',
               'Met': 'CSCC[C@@H](C(=O)O)N',
               'Phe': 'c1ccc(cc1)C[C@@H](C(=O)O)N',
               'Pro': 'C1C[C@H](NC1)C(=O)O',
               'Ser': 'C([C@@H](C(=O)O)N)O',
               'Thr': 'C[C@H]([C@@H](C(=O)O)N)O',
               'Trp': 'c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N',
               'Tyr': 'c1cc(ccc1C[C@@H](C(=O)O)N)O',
               'Arg': 'C(C[C@@H](C(=O)O)N)CNC(=[NH2+])N',
               'Asn': 'C([C@@H](C(=O)O)N)C(=O)N',
               'Gln': 'C(CC(=O)N)[C@@H](C(=O)O)N',
               'Val': 'CC(C)[C@@H](C(=O)O)N',
               'NH2': 'N'}

def extract_table_list(content):
    content_list = [sent for sent in content.split("\n") if len(sent.strip("\n"))>0]
    table_content_list = []
    idx = 0
    while idx < len(content_list):
        sent = content_list[idx]
        if sent.count("|") >= 3:
            temp_range = []
            if idx - 1 >=0 and content_list[idx-1].count("|") <= 2:
                temp_range.append(idx-1)
            temp_range.append(idx)
            idx += 1
            while idx < len(content_list):
                sent = content_list[idx]
                if sent.count("|") >= 3:
                    temp_range.append(idx)
                    idx += 1
                else:
                    break
            if len(temp_range) >= 2:
                temp_range_content = []
                for range_idx in temp_range:
                    temp_range_content.append(content_list[range_idx])
                table_content_list.append(temp_range_content)
                temp_range=[]
        else:
            idx += 1
    return table_content_list

def tabel_check(table_content):
    if len(re.findall(r"[tT]able", table_content[0])) > 0:
        tabel_name = table_content[0]
    else:
        tabel_name = ""

    idx = 1
    col_num = table_content[idx].count("|") - 1
    seq_count = [0 for i in range(col_num)]
    first_line = -1
    first_val = -1
    res = 1
    skip_count = 0
    drop_idx_list = []
    while idx < len(table_content):
        item = table_content[idx]
        if item.count("|") - 1 != col_num:
            drop_idx_list.append(idx)
        if len(re.findall(r"--\s?\|\s?--", item)) > 0:
            first_line = idx - 1
        idx += 1
    if len(drop_idx_list) > len(table_content[first_line:])*0.2:
        res = "表格列无法对齐：{}_{}_{}".format(len(drop_idx_list),len(table_content),first_line)
        return tabel_name, res, -1,drop_idx_list
    elif first_line == -1:
        res = "无法找到首行"
        return tabel_name, res, -1,drop_idx_list
    else:
        return tabel_name, 1, first_line,drop_idx_list
def seq_modify(content):
    content = re.sub(r'<\/?su[bp]>', '', content).strip("\n").strip(" ")
    content = re.sub(r'\[*\(*SEQ.*ID.*\]*\)*','',content)
    content = re.sub(r'^\d+.*\d+$', '', content)
    content_list = content.split('<br>')
    return content_list
def is_seq(content,num_rows):
    content_list = seq_modify(content)
    amino_acids = r"(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)(?:[A-Z\-\(\s]|$)"
    sequence_list = []
    for content in content_list:
        count_AGCTU = sum(1 for char in content if char in {'A', 'C', 'G', 'T', 'U',' '})  # 统计 ACGTU 的数量
        if re.search(r'NH2$|OH$|(?<!form)amide$',content,re.IGNORECASE) :
            sequence_list.append(content)
        # elif set(content).issubset({"A", "C", "G", "T", "U",' ','(',')'}):
        elif content and (count_AGCTU/len(content) > 0.8):  #长字符串识别错误
            pass
        elif len(re.findall(amino_acids, content))>1 and not re.search(r'[,;→]',content) :
            sequence_list.append(content)
        elif set(content).issubset({"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",' ','(',')'}) and ((len(content) > 4 and num_rows >5) or len(content)> 10):
            sequence_list.append(content)
    # if len(sequence_list) > len(content_list)*0.4:
    if sequence_list:
        return True,sequence_list
    else:
        return False,sequence_list
def get_mapping_res(sequence_list):
    sequence_smi =[]
    sequence_three = []
    for sequence in sequence_list:
        print(sequence)
        split_res = sequence_split(sequence)
        match_list = []
        miss_list = []
        for res_item in split_res:
            if res_item in symbol_dict:
                match_list.append({res_item: symbol_dict[res_item]})
            else:
                miss_list.append(res_item)
        if not miss_list and len(split_res) > 1:
            sequence_smi.append(symbol_to_smiles(split_res,symbol_dict))
            sequence_three.append('-'.join(split_res))
    return {
        # "chunk_res": split_res,
        # "match_list": match_list,
        # "miss_list": miss_list,
        'chunk_res':sequence_list,
        'sequence':sequence_three,
        'sequence_smi':sequence_smi
    }

def table2jsonl(table_content,tabel_name,first_line,drop_idx_list):
    new_table_content = []
    for idx in range(len(table_content)):
        if idx not in drop_idx_list:
            new_table_content.append(table_content[idx])
    table_content = new_table_content
    idx = first_line
    col_list = table_content[idx].split("|")
    head_list = col_list[1:-1]
    col_num = len(head_list)
    seq_idx_list,sequence_list = [],[]
    
    seq_count = [0 for i in range(col_num)]
    data_count = [0 for i in range(col_num)]
    # 直接进入第一个有效行
    idx = first_line + 2
    num_rows = len(table_content)
    if idx < num_rows:
        sequence_list = [[None] * len(table_content[idx].split("|")[1:-1]) for _ in range(num_rows)]
    while idx < num_rows:
        item = table_content[idx]
        cell_list = item.split("|")[1:-1]
        for cell_idx,cell_item in enumerate(cell_list):
            # print(cell_item,is_seq(cell_item,num_rows),num_rows,head_list[cell_idx])
            #如果数据所在列的列名包含’name',去除name
            if re.search('name',head_list[cell_idx],re.IGNORECASE):
                sequence = ''
                print(cell_item,is_seq(cell_item,num_rows),num_rows,head_list[cell_idx])
                cell_item = cell_item.split(' ')
                for data in cell_item:
                    status,sequence_new = is_seq(data,10)
                    if status:
                        sequence += ''.join(sequence_new) 
                print(sequence)
                cell_item = sequence
            status,sequence_new = is_seq(cell_item,num_rows)
            if status:
                print(cell_item,is_seq(cell_item,num_rows),num_rows)
                seq_count[cell_idx] += 1
                sequence_list[idx][cell_idx] = sequence_new
            if cell_item.strip() :
                data_count[cell_idx] +=1
        idx += 1
    for seq_idx, item in enumerate(seq_count):
        if item > data_count[seq_idx]*0.8:
            seq_idx_list.append(seq_idx)
    if len(seq_idx_list) == 0:
        return "表格中无法找到序列",[]
     
    idx = first_line + 2
    table_data_list = []
    while idx < num_rows:
        temp_list = []
        for cell_idx, cell_item in enumerate(cell_list):
            key = head_list[cell_idx]
            cell_item = cell_item.strip(" ")
            sequence = sequence_list[idx][cell_idx]
            if cell_idx in seq_idx_list and sequence:
                match_res = get_mapping_res(sequence)
                temp_list.append({
                    key:{"cell_name":cell_item,
                         "match_res":match_res}
                })
            else:
                temp_list.append({
                    key:{"cell_name":cell_item}
                })
        table_data_list.append(temp_list)
        temp_list = []
        idx += 1
    table_data = {
        "table_data_list": table_data_list,
        "seq_idx_list": seq_idx_list,
    }
    return "successful", table_data
fw1 = open("Json_w_seq-312.json","w",encoding="utf-8")
fw2 = open("Json_w_seq-312.txt","w",encoding="utf-8")

fw3 = open("Json_parse_o_seq-312.json","w",encoding="utf-8")
error_table = 0
def parse_single_file(filepath):
    with open(filepath, "r", encoding="utf-8") as fs:
        for item in fs.readlines():
            item = json.loads(item)
            content = item["text"]
            table_list = extract_table_list(content)
            for idx, table_item in enumerate(table_list):
                table_name, check_res, first_line, drop_idx_list = tabel_check(table_item)
                # 表格列无法对齐，或者不存在首行
                if check_res != 1:
                    new_data = copy.deepcopy(item)
                    new_data["table_info"] = {
                        "table_name": table_name,
                        "table_idx": idx,
                        "status": check_res,
                        "drop_idx_list": drop_idx_list,
                        "table_list": "\n".join(table_item),
                        "table_data": {}
                    }
                    # fw2.write(json.dumps(new_data, ensure_ascii=False) + "\n")
                else:
                    status, table_data = table2jsonl(table_item, table_name, first_line, drop_idx_list)
                    new_data = {}
                    if len(table_data) == 0:
                        new_data = copy.deepcopy(item)
                        new_data["table_info"] = {
                            "table_name": table_name,
                            "table_idx": idx,
                            "status": status,
                            "drop_idx_list": drop_idx_list,
                            "table_list": "\n".join(table_item),
                            "table_data": {}
                        }
                        fw3.write('page: '+str(item['page'])+'\n')
                        fw3.write('file: '+str(item['file'])+'\n')
                        fw3.write(new_data["text"]+ "\n")
                    else:
                        new_data = copy.deepcopy(item)
                        new_data["table_info"] = {
                            "table_name": table_name,
                            "table_idx": idx,
                            "status": status,
                            "drop_idx_list": drop_idx_list,
                            "table_list": "\n".join(table_item),
                            "table_data": table_data
                        }
                        fw1.write(json.dumps(new_data, ensure_ascii=False) + "\n")
                        fw2.write('page: '+str(item['page'])+'\n')
                        fw2.write('file: '+str(item['file'])+'\n')
                        fw2.write(item['text']+ "\n")
                        

# base_dir = "../final_res/"
# for file in os.listdir(base_dir):
#     file = "{}/{}".format(base_dir, file)
#     parse_single_file(file)

# file = '../final_res/step4_tabel2Json_w_seq.jsonl'

file = '../final_res/step4_tabel2Json_parse_o_seq.jsonl'
parse_single_file(file)
