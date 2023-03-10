import copy
import os

from Bio.PDB import PDBIO, Superimposer

from RNA_normalizer.msgs import show
from RNA_normalizer.utils import erf

BIN_DIR = os.getcwd()


class PDBComparer:
    BACKBONE_ATOMS = [
        "C1'",
        "C2'",
        "C3'",
        "C4'",
        "C5'",
        "O2'",
        "O3'",
        "O4'",
        "O5'",
        "OP1",
        "OP2",
        "P",
    ]
    HEAVY_ATOMS = [
        "C2",
        "C4",
        "C5",
        "C6",
        "C8",
        "N1",
        "N2",
        "N3",
        "N4",
        "N6",
        "N7",
        "N9",
        "O2",
        "O4",
        "O6",
    ]
    ALL_ATOMS = BACKBONE_ATOMS + HEAVY_ATOMS

    RMSDD_ATOMS = ["C4", "C8", "P", "C1'"]

    def __init__(self):
        pass

    def mcq(self, f1, f2):
        cmd = (
            "java -cp %s/mcq.ws.client-0.0.1-SNAPSHOT-jar-with-dependencies.jar pl.poznan.put.mcq.ws.client.Global -m %s -t %s >mcq.log"
            % (BIN_DIR, f1, f2)
        )
        os.system(cmd)
        try:
            v = float(open("mcq.log").read().strip())
        except Exception:
            v = 0
        return v

    def gdt(self, f1, f2):
        cmd = "java -jar %s/gdt.jar %s %s >gdt.log" % (BIN_DIR, f2, f1)
        # ~ print cmd
        os.system(cmd)
        try:
            x = open("gdt.log").read().strip().split("\n")[1].split(",")[-1]
            if x == "NaN":
                return 0
            v = float(x)
        except Exception:
            v = 0
        return v

    def rmsd(self, src_struct, trg_struct, fit_pdb=None):
        # for each model in the reference
        src_residues = src_struct.res_sequence()
        trg_residues = trg_struct.res_sequence()

        atoms = self._get_atoms_struct(
            PDBComparer.ALL_ATOMS, src_residues, trg_residues
        )

        if not atoms is None:
            (src_atoms, trg_atoms) = atoms
        else:
            return None

        # compute the rmsd value and apply it to the target structure
        sup = Superimposer()
        sup.set_atoms(src_atoms, trg_atoms)

        # we copy the fit_struct to leave the target struct unmodified for posterior processing
        fit_struct = copy.deepcopy(trg_struct.struct)
        sup.apply(fit_struct.get_atoms())

        # save the fitted structure
        if not fit_pdb is None:
            io = PDBIO()
            io.set_structure(fit_struct)
            io.save(fit_pdb)

        return sup.rms

    # From Hajdin et al., RNA (7) 16, 2010
    def pvalue(self, m, N, param):
        if param == "+":
            a = 5.1
            b = 15.8
        elif param == "-":
            a = 6.4
            b = 12.7
        else:
            show("FATAL", "Wrong p-value parameter '%s'. Expected '+' or '-'" % param)

        RMSD = a * (N**0.41) - b

        Z = (m - RMSD) / 1.8

        pv = (1.0 + erf(Z / (2**0.5))) / 2.0

        return pv

    def INF(self, src_struct, trg_struct, type):
        (P, TP, FP, FN) = (0, 0, 0, 0)

        for stype, sb1, sb2, sextra in src_struct.get_interactions(type):
            P += 1
            found = False
            for ttype, tb1, tb2, textra in trg_struct.get_interactions(type):
                if (
                    (stype == ttype)
                    and (sb1 == tb1)
                    and (sb2 == tb2)
                    and (sextra == textra)
                ):
                    found = True
                    break

            if found:
                # print "TP>", (stype, sb1, sb2, sextra)
                TP += 1
            else:
                # print "FN>", (stype, sb1, sb2, sextra)
                FN += 1

        for ttype, tb1, tb2, textra in trg_struct.get_interactions(type):
            found = False
            for stype, sb1, sb2, sextra in src_struct.get_interactions(type):
                if (
                    (stype == ttype)
                    and (sb1 == tb1)
                    and (sb2 == tb2)
                    and (sextra == textra)
                ):
                    found = True
                    break

            if not found:
                FP += 1
            # print "FP>", (ttype, tb1, tb2, textra)

        if TP == 0 and (FP == 0 or FN == 0):
            INF = -1.0
        else:
            PPV = float(TP) / (float(TP) + float(FP))
            STY = float(TP) / (float(TP) + float(FN))
            INF = (PPV * STY) ** 0.5

        # print "##>", INF, P, TP, FP, FN
        return INF

    def DP(self, src_struct, trg_struct, template_txt, dname, dp_script):
        # prepare the config file
        txt = ""
        txt += "matrix=True\n"
        txt += "quiet_err = True\n"
        txt += "out_dir = '%s'\n" % dname
        txt += "ref_model = ('%s', 0)\n" % src_struct.pdb_file
        txt += "cmp_model = [('%s', 0)]\n" % trg_struct.pdb_file

        aligns = self._build_dp_alignments(src_struct, trg_struct)
        aligns_txt = []

        for align in aligns:
            aligns_txt.append(
                "('%s', %s, '%s', %s, %s)"
                % (align[0], align[1], align[2], align[3], align[4])
            )

        txt += "aligns = [%s]\n" % (", ".join(aligns_txt))
        txt += template_txt

        fname_cfg = "%s.cfg" % (trg_struct.pdb_file)
        fname_log = "%s.log" % (trg_struct.pdb_file)
        open(fname_cfg, "w").write(txt)

        # runs the DP generator
        os.system("python %s -c %s > %s" % (dp_script, fname_cfg, fname_log))

    def VARNA(self, src_struct, trg_struct, algorithm="radiate"):
        edges = {"W": "wc", "S": "s", "H": "h"}

        data = {}
        data["sequenceDBN"] = src_struct.raw_sequence()
        data["structureDBN"] = "." * len(src_struct.raw_sequence())

        aux_bps = []

        for stype, sb1, sb2, sextra in src_struct.get_interactions("PAIR"):
            color = "#FF0000"
            for ttype, tb1, tb2, textra in trg_struct.get_interactions("PAIR"):
                if (
                    (stype == ttype)
                    and (sb1 == tb1)
                    and (sb2 == tb2)
                    and (sextra == textra)
                ):
                    color = "#00FF00"
                    break
            aux_bps.append(
                "(%d,%d):color=%s,edge5=%s,edge3=%s,stericity=%s"
                % (
                    sb1 + 1,
                    sb2 + 1,
                    color,
                    edges[sextra[0]],
                    edges[sextra[1]],
                    sextra[2:],
                )
            )

        data["auxBPs"] = ";".join(aux_bps)
        data["algorithm"] = algorithm

        return data

    def _get_atoms_residue(self, atom_list, src_res, trg_res):
        src_atom_list = []
        trg_atom_list = []

        src_atom_list_tmp = list(filter(lambda a: a.get_name() in atom_list, src_res))
        trg_atom_list_tmp = list(filter(lambda a: a.get_name() in atom_list, trg_res))

        # for each atom in reference
        for src_atom in src_atom_list_tmp:
            found = False
            src_name = src_atom.get_full_id()[4][0]

            # search for an atom with the same name in the comparison
            for trg_atom in trg_atom_list_tmp:
                trg_name = trg_atom.get_full_id()[4][0]

                # if the atom was found keep it and jump to the next
                if src_name == trg_name:
                    src_atom_list.append(src_atom)
                    trg_atom_list.append(trg_atom)
                    found = True
                    break

            if not found:
                show(
                    "WARNING",
                    "Atom %s from residue %s not found in target atom list"
                    % (src_name, src_res.id),
                )
        return (src_atom_list, trg_atom_list)

    def _get_atoms_struct(self, atom_list, src_residues, trg_residues):
        src_atoms = []
        trg_atoms = []

        if len(src_residues) != len(trg_residues):
            show("ERROR", "Different number of residues!")
            return None

        for src_res, trg_res in zip(src_residues, trg_residues):
            (sa, ta) = self._get_atoms_residue(atom_list, src_res, trg_res)

            src_atoms.extend(sa)
            trg_atoms.extend(ta)
        # 'print('%d %d'%(len(src_atoms),len(trg_atoms)))
        return (src_atoms, trg_atoms)

    def _build_dp_alignments(self, src_struct, trg_struct):
        aligns = []

        (schain, tchain) = ("", "")
        (spos, tpos) = (-1, -1)
        count = 0
        item = None

        for i, j in zip(src_struct.res_seq, trg_struct.res_seq):
            (sres, tres) = (src_struct.res_list[i], trg_struct.res_list[j])

            # if the numbering or the chain change
            if (
                (sres.pos != (spos + 1))
                or (tres.pos != (tpos + 1))
                or (sres.chain != schain)
                or (tres.chain != tchain)
            ):
                if count > 0:
                    item[4] = count
                    aligns.append(item)

                (schain, tchain) = (sres.chain, tres.chain)
                item = [sres.chain, sres.pos, tres.chain, tres.pos, None]
                count = 1
            else:
                count += 1

            (spos, tpos) = (sres.pos, tres.pos)

        # adds the last alignment
        if count > 0:
            item[4] = count
            aligns.append(item)

        return aligns
