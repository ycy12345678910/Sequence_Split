
from rdkit import Chem
from rdkit.Chem import AllChem as Chem

def symbol_to_smiles(symbol_split,symbol_dict):
    carboxyl_smiles = symbol_dict[symbol_split[0]]
    aa_num = 0
    for symbol in symbol_split[1:]:
        carboxyl_mol = Chem.MolFromSmiles(carboxyl_smiles)  # 创建前一个氨基酸
        carboxyl_count = carboxyl_smiles.count("=O") + carboxyl_smiles.count("O=") -aa_num
        amine_smiles = symbol_dict[symbol] 
        amine_mol = Chem.MolFromSmiles(amine_smiles)  # 创建下一个氨基酸
        amine_count = amine_smiles.count('N') 
        if symbol == 'Pim':
            reaction_smarts = "[C:0][C:1](=[O:2])[O:3].[c:4]([n:5])[#7:6]>>[C:0][c:4]([n:5])[#7:6].[C:1](=[O:2])[O:3]"
        elif amine_count == carboxyl_count == 1:
            reaction_smarts = "[C:1](=[O:2])[O:3].[N:4]>>[C:1](=[O:2])-[N:4].[O:3]"  
        elif amine_count ==1:
            reaction_smarts = "[N:1][C:2][C:3](=[O:4])[O:5].[N:6]>>[N:1][C:2][C:3](=[O:4])-[N:6].[O:5]"  
        elif carboxyl_count ==1:
            reaction_smarts = "[C:1](=[O:2])[O:3].[N:4][C:5][C:6](=[O:7])[O:8]>>[C:1](=[O:2])-[N:4][C:5][C:6](=[O:7])[O:8].[O:3]" 
        else:
            reaction_smarts = "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1][C:2][C:3](=[O:4])-[N:6][C:7][C:8](=[O:9])[O:10].[O:5]"  
        
        reaction = Chem.ReactionFromSmarts(reaction_smarts)  # 创建反应
        # 3. 执行反应，连接两个分子
        products = reaction.RunReactants([carboxyl_mol,amine_mol])
        aa_num += 1
        if products:
            product_smiles = Chem.MolToSmiles(products[0][0])
            carboxyl_smiles = product_smiles
        else:
            return '*'*10,symbol,amine_count,carboxyl_count
    return product_smiles