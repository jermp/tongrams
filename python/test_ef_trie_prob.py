import tongram
import time

model_handler1 = tongram.load_ef_trie_prob_backoff("sample_text_arpa_ef_trie.prob_backoff.bin")

start = time.time()
score = tongram.score_ef_trie_prob_backoff(model_handler1, "Take this computer")
end = time.time()
print 'Score: '+str(score)
print(str(1000*(end - start))+' ms')