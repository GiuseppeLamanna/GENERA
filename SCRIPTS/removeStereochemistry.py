import sys

fh = open(sys.argv[1], "r")
chromo = fh.read().split("\n")[0]
fh.close()


out = chromo.replace("/", "").replace("\\", "").replace("@", "") +"\n"


fh = open(sys.argv[2], "w")
fh.write(out)
fh.close()
