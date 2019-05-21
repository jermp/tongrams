import tongram

model_handler1 = tongram.load("../debug/pef_trie.count.out_1")
model_handler2 = tongram.load("../debug/pef_trie.count.out_2")
print 'done loading model'

print 'Lookup: '+str(tongram.lookup(model_handler1, "want this"))
print 'Lookup: '+str(tongram.lookup(model_handler2, "want this"))
