
maxthreads = 384
numthreads = 384
numblocks = 1024
numiters = 16

defnumevts = numthreads * numblocks * numiters
curnumevts = defnumevts


# while (numthreads > 0 and defnumevts == curnumevts):
#     print numblocks, numthreads, numiters
#     numthreads /= 2
#     numblocks *= 2
#     curnumevts = numblocks * numthreads * numiters
#
# numthreads = 1
# numblocks = defnumevts / numiters
# curnumevts = numblocks * numthreads * numiters
#
#
# while (numthreads <= maxthreads and defnumevts == curnumevts):
#     print numblocks, numthreads, numiters
#     numthreads *= 2
#     numblocks /= 2
#     curnumevts = numblocks * numthreads * numiters

cnt = 0


blocklist = [2**x for x in range(0, 12)]
threadlist = [1, 3, 6, 12, 24, 48, 96, 192, 384]
maxiter = 8
numevts = maxiter * blocklist[-1] * threadlist[-1]
print numevts
print blocklist
for block in blocklist:
    for thread in threadlist:
        iter = numevts / block / thread
        print block, thread, iter
