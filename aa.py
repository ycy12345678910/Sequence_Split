import re

def remove_hyphens(s):
    return re.sub(r'(?<=[a-z])-(?=[a-z])', '', s)

s = "(2-Methyl-3-(3-indolyl)propionyl)-Lys(ε-N-(2-methylphenyl)amino-carbonyl)-Asp-(NMe)PheNH₂"
result = remove_hyphens(s)
print(result)
