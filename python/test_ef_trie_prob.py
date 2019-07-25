import tongram
import time

model_handler1 = tongram.load("sample_text_arpa_ef_trie.prob_backoff.bin")

start = time.time()
score = tongram.score(model_handler1, "Take this computer")
end = time.time()
print 'Score: '+str(score)
print(str(1000*(end - start))+' ms')
