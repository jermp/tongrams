import tongram

tongram.load("../data_pef_trie.count.out")
print 'done loading model'

print 'Lookup: '+str(tongram.lookup("this is a"))

