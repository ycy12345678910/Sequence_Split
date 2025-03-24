from rdkit import Chem
from rdkit.Chem import AllChem as Chem
import re
import pandas as pd

symbol_smiles = {'Asp':'C([C@@H](C(=O)O)N)C(=O)O','Cys':'C([C@@H](C(=O)O)N)S','Dap':'C(C(C(=O)O)N)N','Glu':'C(CC(=O)O)[C@@H](C(=O)O)N','Ser':'C([C@@H](C(=O)O)N)O',
                 'Lys':'C(CCN)C[C@@H](C(=O)O)N','Gly':'C(C(=O)O)N','Trp':'c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N'}
def name2smiles(side_modify):
    if side_modify in yl_smiles.keys():
        modify_smiles = yl_smiles[side_modify]
    else: 
        side_modify_new = '('+ side_modify +')chloride-37Cl' #对取代基后添加[37Cl]定位
        #通过ChemDraw转化为smiles  C:\Python32\python.exe .\name2smi.py "side_modify_new"
        modify_smiles = input()
        yl_smiles[side_modify] = modify_smiles
    return modify_smiles
    
def reaction_sidechain(AA_base,modify_smiles):
    AA_mol = Chem.MolFromSmiles(symbol_smiles[AA_base])
    modify_mol = Chem.MolFromSmiles(modify_smiles)
    if AA_base in ['Lys','Ser','Dap','Cys']:  #氨基酸直接接取代基
        reaction_smarts = "[C;H2:0][N,O,S:1].[37Cl:2][A,a:3] >> [C;H2:0][N,O,S:1][A,a:3].[37Cl:2]"
    elif AA_base in ['Asp','Glu']:  #脱水缩合
        reaction_smarts = "[C;H2:0][C,H0:4][O;H1:1].[37Cl:2][A,a:3] >> [C;H2:0][C,H0:4][A,a:3].[37Cl:2].[O;H1:1]" 
    elif AA_base == 'Gly':
        reaction_smarts = "[C;H2:0][N;H2:1].[37Cl:2][A,a:3] >> [N;H2:1][C;H1:0][A,a:3].[37Cl:2]"
    else:
        print("请输入反应类型,eg:[C;H2:0][N,O,S:1].[37Cl:2][A,a:3] >> [C;H2:0][N,O,S:1][A,a:3].[37Cl:2]")
        reaction_smarts = input()
    reaction = Chem.ReactionFromSmarts(reaction_smarts)  # 创建反应
    AA_mol = reaction.RunReactants([AA_mol,modify_mol])
    return AA_mol
def modify_sidechain(matches_1):
    AA_base = matches_1.group(1)
    side_modify = matches_1.group(2)
    if side_modify.startswith(("NH-", "O-","S-")):  #判断是否是NH-/O-/S-开头
        side_modify_new = re.sub(r"^(NH-|O-|S-)",'',side_modify)
        modify_smiles = name2smiles(side_modify_new)
        modify_smiles = modify_smiles.replace("[37Cl]", f'{side_modify[0]}[37Cl]') #在([37Cl])前面插入修饰元素
    elif side_modify.endswith('yl'):
        modify_smiles = name2smiles(side_modify)
    else:
        print(side_modify,'其它情形,请输入smiles,以[37Cl]定义连接位点,eg:CCCCCCCO[37Cl]')
        modify_smiles = name2smiles(side_modify)
    AA_mol = reaction_sidechain(AA_base,modify_smiles) 
    if AA_mol:
        AA_smiles = Chem.MolToSmiles(AA_mol[0][0])
        return AA_smiles
def reaction_mainchain(AA_base,modify):
    if modify in ['NMe','N-methyl']:
        reaction_smarts = "[N;!H0:0][C:1][C;H0:2][O;H1:3].[37Cl:4][A,a:5] >>[A,a:5][N;!H0:0][C:1][C;H0:2][O;H1:3].[37Cl:4]"
        modify_smiles = 'C[37Cl]'
    elif modify.startswith("alpha"):
        modify = modify.replace('alpha-','')
        reaction_smarts = "[N;!H0:0][C;!H0:1][C;H0:2][O;H1:3].[37Cl:4][A,a:5] >>[N;!H0:0][C:1]([A,a:5])[C;H0:2][O;H1:3].[37Cl:4]"
        modify_smiles = yl_smiles[modify]
    else:
        print(modify,'其它情形,请输入smiles,以[37Cl]定义连接位点,eg:CCCCCCCO[37Cl]')
        modify_smiles = name2smiles(modify)
    AA_mol = Chem.MolFromSmiles(symbol_smiles[AA_base])
    modify_mol = Chem.MolFromSmiles(modify_smiles)
    reaction = Chem.ReactionFromSmarts(reaction_smarts)  # 创建反应
    AA_mol = reaction.RunReactants([AA_mol,modify_mol])
    return AA_mol
def modify_mainchain(matches_2):
    print(matches_2.group(1))
    AA_base = re.sub(r"^(-)",'',matches_2.group(2))
    modify = matches_2.group(1)
    AA_mol = reaction_mainchain(AA_base,modify)
    if AA_mol:
        AA_smiles = Chem.MolToSmiles(AA_mol[0][0])
        return AA_smiles

AA_list = ['Asp(1-heptanol)','Asp(NH-hexyl)','Cys(S-heptyl)','Dap(octanoyl)','Glu(O-hexyl)','Ser(n-octanoyl)','Lys(biotinyl)','Gly(myristyl)','(alpha-Methyl)-Trp','(NMe)Asp',
           'Lys(2-(2-(2-(4-(Hexadecanoylamino)-4(S)-carboxybutyrylamino)ethoxy)ethoxy)acetyl)','t-BOC-Trp']
df_symbol = pd.read_csv('yl_smiles.csv') 
yl_smiles = dict(zip(df_symbol['yl'],df_symbol['SMILES']))
for AA in AA_list:
    print(AA)
    matches_1 = re.search(r'([0-9a-zA-Z\-]+)\((.*)\)',AA)  #侧链修饰
    matches_2 = re.search(r'\((.*)\)([D-]*[0-9a-zA-Z]{3})|(.+)\-([D-]*[0-9a-zA-Z]{3})',AA)  #主链修饰
    # matches_2 = re.search(r'\((.*)\)([0-9a-zA-Z\-]+)',AA)  #主链修饰
    print(matches_1,matches_2)
    if matches_1:
        AA_smiles = modify_sidechain(matches_1)
        print(AA_smiles)
    elif matches_2:
        AA_smiles = modify_mainchain(matches_2)
        print(AA_smiles,'2'*20)
df = pd.DataFrame(list(yl_smiles.items()), columns=['yl', 'SMILES'])
df.to_csv("yl_smiles.csv", index=False,)  
