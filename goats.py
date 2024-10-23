import mygene as mg

mg = mg.MyGeneInfo()

out = mg.querymany(['GAPDH', 'CCDC83'], scopes='symbol', species='human', fields=['go'])

print(out)
