import sourcerandom

#RAND_GEN = sourcerandom.SourceRandom(source=sourcerandom.OnlineRandomnessSource.RANDOM_ORG)
RAND_GEN = sourcerandom.SourceRandom(source=sourcerandom.OnlineRandomnessSource.QRNG_ANU) 
# Initialize the generator - use the same generator to take advantage of caching

x = pow(2, 512)

f = open('random_number.txt', 'w')
print((hex(RAND_GEN.randint(0, x))).upper(), file=f) 
    # Every method from random package is included