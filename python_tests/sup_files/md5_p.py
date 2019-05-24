import hashlib, sys


def mdd5(filename):
    hasher = hashlib.md5()
    f = open(filename, "rb")
    content = f.read()
    hasher.update(content)
    #print(hasher.hexdigest())
    f.close()
    return hasher.hexdigest()
