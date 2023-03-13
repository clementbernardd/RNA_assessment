import math
import sys

from Bio.PDB import PDBParser, Vector

from lib.rna_assessment.RNA_normalizer.mcannotate import MCAnnotate
from lib.rna_assessment.RNA_normalizer.msgs import show
from lib.rna_assessment.RNA_normalizer.structures.residue import Residue


class PDBStruct(object):
    def __init__(self, mc_annotate_bin):
        self._pdb_file = None
        self._struct = None
        self._res_list = []
        self._res_seq = []
        self._res_index = {}
        self._interactions = []
        self.mc_annotate_bin = mc_annotate_bin

    # self._brackets = []
    # self._wcpairs = []

    def load(self, pdb_file, index_name=None):
        self._pdb_file = pdb_file

        ok = self._load_struct()

        if ok and index_name is not None:
            ok = self._load_index(index_name)
        else:
            ok = self._load_index2()

        if ok:
            ok = self._load_annotations_3D()
        # ~ print pdb_file,self._interactions,index_name
        # if( ok and not fbrackets is None ):
        # 	ok = self._load_brackets( fbrackets )
        return ok

    def raw_sequence(self):
        seq = ""
        for ndx in self._res_seq:
            seq += self._res_list[ndx].nt

        return seq

    def res_sequence(self):
        result = []
        for ndx in self._res_seq:
            result.append(self._res_list[ndx].res)

        return result

    def get_interactions(self, type="ALL"):
        if type == "ALL":
            # "ALL": returns all interactions
            return self._interactions
        elif type in ("PAIR"):
            # "PAIR": returns all pairs irrespective of their type
            return list(
                filter(lambda x: x[0] in ("PAIR_2D", "PAIR_3D"), self._interactions)
            )
        elif type in ("PAIR_2D", "PAIR_3D", "STACK"):
            # "PAIR_2D", "PAIR_3D", "STAK": returns the interactions of the specified type
            return list(filter(lambda x: x[0] == type, self._interactions))
        else:
            show(
                "FATAL",
                "Wrong interaction type '%s' expected: 'ALL', 'PAIR', 'PAIR_2D', 'PAIR_3D' or 'STACK'"
                % type,
            )

    # --- properties ---
    def struct_get(self):
        return self._struct

    def res_seq_get(self):
        return self._res_seq

    def res_list_get(self):
        return self._res_list

    def pdb_file_get(self):
        return self._pdb_file

    def rad_gir(self):
        rmean = Vector([0.0, 0.0, 0.0])
        count = 0
        for res in self._res_list:
            for a in res.res:
                rmean += a.get_vector()
                count += 1

        rmean = rmean / float(count)

        rsum = 0.0
        for res in self._res_list:
            for a in res.res:
                rsum += (a.get_vector() - rmean) * (a.get_vector() - rmean)

        return math.sqrt(rsum / count)

    # def brackets_get(self):
    # 	return self._brackets

    struct = property(struct_get)
    res_seq = property(res_seq_get)
    res_list = property(res_list_get)
    pdb_file = property(pdb_file_get)

    # brackets = property( brackets_get )
    # ---

    def _load_struct(self):
        parser = PDBParser(QUIET=True)
        self._struct = parser.get_structure("struct", self._pdb_file)

        if len(self._struct) > 1:
            show(
                "WARNING",
                "%d models found. Only the first will be used!" % (len(self._struct)),
            )

        self._res_list = []
        self._res_seq = []
        self._res_index = {}

        # gets only the first model
        model = self._struct[0]
        count = 0
        for chain in model.child_list:
            for res in chain.child_list:
                new_residue = Residue(chain.id, res.id[1], res.resname.strip(), res)

                self._res_list.append(new_residue)
                self._res_seq.append(count)
                self._res_index[new_residue.key()] = [count, None]

                count += 1

        return True

    def _load_index(self, index_name):
        self._res_seq = []
        entries = []
        for row in open(index_name).read().split("\n"):
            row = row.strip()
            if (not row.startswith("#")) and (row != ""):
                entries.extend(map(lambda row: row.split(":"), row.split(",")))

        for entry in entries:
            if len(entry) != 3:
                show("ERROR", "Bad index entry: '%s'" % entry)
                return False

            chain = entry[0]
            pos = int(entry[1])
            count = int(entry[2])

            # get the index position
            ndx = self._get_index(chain, pos, 0)

            if ndx is None:
                return False

            # get the positions
            for i in range(ndx, ndx + count):
                if i >= len(self._res_list):
                    show("ERROR", "Bad count %d in index entry: '%s'" % (count, entry))
                    return False

                if self._res_list[i].chain != chain:
                    show(
                        "ERROR",
                        "Position %d in index entry: '%s' is outside the chain"
                        % (i, entry),
                    )
                    return False
                self._res_seq.append(i)

                # update the index with the rank of the residue
                self._res_index[self._res_list[i].key()][1] = len(self._res_seq) - 1
        return True

    def _load_index2(self):
        self._res_seq = []
        for i in range(0, len(self._res_list)):
            self._res_seq.append(i)
            self._res_index[self._res_list[i].key()][1] = len(self._res_seq) - 1
        return True

    def _load_annotations_3D(self):
        self._interactions = []
        mca = MCAnnotate(self.mc_annotate_bin)
        mca.load(self._pdb_file)
        # ~ print mca.interactions
        for (
            type,
            chain_a,
            pos_a,
            nt_a,
            chain_b,
            pos_b,
            nt_b,
            extra1,
            extra2,
            extra3,
        ) in mca.interactions:
            # get the rank of the first position of the pair
            rank_a = self._get_index(chain_a, pos_a, 1)
            rank_b = self._get_index(chain_b, pos_b, 1)

            if (rank_a is None) or (rank_b is None):
                continue
            # ~ return False

            if type == "STACK":
                extra = extra1
            else:
                extra = "%s%s" % (extra1, extra2)
            self._interactions.append(
                (type, min(rank_a, rank_b), max(rank_a, rank_b), extra)
            )

    def _get_index(self, chain, pos, field):
        key = "%s:%s" % (chain, pos)
        data = self._res_index.get(key, None)[field]
        if data is None and field == 0:
            sys.stderr.write("ERROR	Bad index key: '%s'\n" % key)

        return data
