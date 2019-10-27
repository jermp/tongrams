import tongram

model = tongram.load("../build/pef_trie.count.bin")
print 'done loading model'
print 'Lookup: '+str(tongram.lookup(model, "this is a"))
