from Bio import trie
trieobj = trie.trie()
trieobj["hello"] = 5
trieobj["hello"] += 1
trieobj["abc"] = dict()
trieobj["abc"]["zzz"] = 5
print trieobj["hello"]
print trieobj.keys()
print trieobj.values()
#there is not .items() method

