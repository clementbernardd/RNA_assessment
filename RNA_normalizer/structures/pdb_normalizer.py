from lib.rna_assessment.RNA_normalizer.msgs import show


class PDBNormalizer:
    MAX_ERRORS = 5

    def __init__(self, fres_list, fatoms_list):
        self._load_res_list(fres_list)
        self._load_atom_list(fatoms_list)

    def parse(self, finput, foutput):
        # state variables for the parse process
        self._in_model = False
        self._in_atom = False

        self._chain_found = False
        self._row_count = 0
        self._ok = True

        out_txt = ""

        fi = open(finput)

        for row in fi:
            self._row_count += 1

            row = row.strip()

            rec_name = row[:6]

            if rec_name == "MODEL ":
                row = self.parse_model(row)
                row = ""
            elif rec_name == "ENDMDL":
                row = self.parse_endmdl(row)
                row = ""
            elif rec_name[:3] == "TER":
                row = self.parse_ter(row)
            elif rec_name in ("ATOM  ", "HETATM"):
                row = self.parse_atom(row)
            else:
                continue

            if row != "":
                out_txt += row + "\n"

        fi.close()

        if self._in_atom:
            out_txt += "TER\n"

        if self._ok:
            open(foutput, "w").write(out_txt)

        return self._ok

    def parse_model(self, row):
        if self._in_model:
            self.show_err("'ENDMDL' not found.")

        if self._in_atom:
            self.show_err("Missing 'MODEL' before 'ATOM' declaration.")

        self._in_model = True
        return row

    def parse_endmdl(self, row):
        if not self._in_model:
            show("Warning", "Missing 'MODEL' declaration.")
        # ~ self.show_err( "Missing 'MODEL' declaration." )

        if self._in_atom:
            show("Warning", "Missing 'TER' declaration.")
        # ~ self.show_err( "Missing 'TER' declaration." )
        self._in_model = False
        self._in_atom = False
        return "ENDMDL"

    def parse_ter(self, row):
        result = ""

        if self._in_atom:
            result = "TER"

        self._in_atom = False
        return result

    def parse_atom(self, row):
        # get all the fields from the line
        serial = row[6:11]
        name = row[12:16].strip()
        altLoc = row[16]
        resName = row[17:20].strip()
        chainID = row[21]
        resSeq = int(row[22:26])
        iCode = row[26]
        x = row[30:38]
        y = row[38:46]
        z = row[46:54]
        occupancy = row[54:60]
        tempFactor = row[60:66]
        element = row[76:78]
        charge = row[78:80]

        # check residue name
        resName_norm = self._res_list.get(resName, None)
        if resName_norm is None:
            self.show_err("Unknown residue name: '%s'." % resName)
            return ""
        elif resName_norm == "-":
            return ""

        resName = resName_norm

        # check atom name
        name_norm = self._atom_list.get(name, None)
        if name_norm is None:
            self.show_err("Unknown atom name: '%s' in residue'%s'" % (name, resName))
            return ""
        elif name_norm == "-":
            return ""

        name = name_norm.ljust(3)

        # check chainID
        if chainID == " ":
            if self._chain_found:
                self.show_err("One of the chains is missing!")
            else:
                chainID = "A"
        else:
            self._chain_found = True

        # check occupancy
        # if( occupancy == "" ):

        # This is the only that works with MolProbity
        occupancy = "  1.00"

        # check tempFactor
        if tempFactor == "":
            tempFactor = "  0.00"

        # check element
        if element == "":
            element = name[0]

        self._in_atom = True
        return "ATOM  %5s  %3s%s%3s %s%4d%s   %8s%8s%8s%6s%6s		  %2s%2s" % (
            serial,
            name,
            altLoc,
            resName,
            chainID,
            resSeq,
            iCode,
            x,
            y,
            z,
            occupancy,
            tempFactor,
            element,
            charge,
        )

    def show_err(self, msg):
        show("ERROR", "Line %d: %s\n" % (self._row_count, msg))
        self._ok = False

    def _load_res_list(self, fres_list):
        # read residues list
        pairs = map(
            lambda x: x.split(),
            filter(
                lambda row: not row.startswith("#"),
                open(fres_list).read().strip().split("\n"),
            ),
        )

        self._res_list = {}
        for name, nt in pairs:
            self._res_list[name] = nt

    def _load_atom_list(self, fatoms_list):
        # read residues list
        pairs = map(
            lambda x: x.split(),
            filter(
                lambda row: not row.startswith("#"),
                open(fatoms_list).read().strip().split("\n"),
            ),
        )

        self._atom_list = {}
        for name, name_norm in pairs:
            self._atom_list[name] = name_norm
