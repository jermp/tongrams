import tongrams
import time

count_model_filename = "../build/ef_trie.count.bin"
count_model = tongrams.CountModel(count_model_filename)
ngrams = ["This is a", "function call", "compiler", "this is a"]
for ngram in ngrams:
    print("count of '" + ngram + "' is: " + str(count_model.lookup(ngram)))

def perplexity(log10_prob, words):
    return pow(10, -log10_prob / words)

prob_model_filename = "../build/ef_trie.prob_backoff.bin"
prob_model = tongrams.ProbModel(prob_model_filename)

courpus_filename = "../test_data/sample_text"
start = time.time()
score = prob_model.score_corpus(courpus_filename)
end = time.time()
print("Score: " + str(score))
print(str((end - start)) + " sec")
print("Perplexity is: " + str(perplexity(score[0], score[1])))

sentence = "Having a little flexibility on that issue would go a long way to putting together a final package ."
start = time.time()
score = prob_model.score_sentence(sentence, len(sentence))
end = time.time()
print("Score: " + str(score))
print(str((end - start)) + " sec")
print("Perplexity is: " + str(perplexity(score[0], score[1])))