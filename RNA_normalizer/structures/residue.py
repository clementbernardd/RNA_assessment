#
# get the sequence list from a pdb either raw or indexed
#
class Residue:
    def __init__(self, chain, pos, nt, res):
        self.chain = chain
        self.pos = pos
        self.nt = nt
        self.res = res

    def key(self):
        return "%s:%s" % (self.chain, self.pos)

    def __str__(self):
        return "%s:%s:%s > %s" % (self.chain, self.pos, self.nt, self.res)
