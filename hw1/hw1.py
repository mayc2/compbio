from Bio import SeqIO

for n in SeqIO.parse("sequence.fasta","fasta"):
	next

#test
print(n.seq)

a = n.seq.count("A")
c = n.seq.count("C")
t = n.seq.count("T")
g = n.seq.count("G")
at = a + t
cg = c + g

print(a,c,t,g, at,cg)

l1 = 500
l2 = 1000
l3 = len(n.seq)/20
print len(n.seq)
print(l1,l2,l3)

