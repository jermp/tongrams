import tongram
import time

model = tongram.load("../build/pef_trie.prob_backoff.bin")
start = time.time()
score = tongram.score(model, "../test_data/sample_text.LESSER")
end = time.time()
print 'Score: '+str(score)
print(str(1000*(end - start))+' ms')