##
# File:    SwissModelPdbToCifConverter.py
# Author:  Dennis Piehl
# Date:    01-Oct-2021
#
# Update:
#
#
# References:
#  - Code originally based on: https://github.com/rcsb/modbase_utils
##
"""
Processors for converting SWISS-MODEL PDB files to mmCIF format.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import datetime
import logging

from rcsb.utils.modbase_utils.modbase_pdb_to_cif import CifLoop, three_to_one, one_to_three
from rcsb.utils.io.FileUtil import FileUtil

logger = logging.getLogger(__name__)


class SwissModelPdbToCifConverter:
    """Class for reading SWISS-MODEL PDB files and writing out mmCIF files."""

    def __init__(self, **kwargs):
        self.__cachePath = kwargs.get("cachePath", "./CACHE-insilico3d-models")
        self.__fU = FileUtil(workPath=self.__cachePath)
        #
        # Initialize all attributes
        self.pdbRead = False
        self.fh = None
        self.organism = kwargs.get("organism", None)
        self.uniProtId = kwargs.get("uniProtId", None)
        self.title, self.expdta, self.author, self.revDate = None, None, None, None
        self.jrnl, self.remarks, self.atoms = [], [], []
        self.align = False
        self.todaysDate = datetime.datetime.today().strftime("%Y-%m-%d")
        self.modelFileName, self.modelNameBase = None, None
        self.atomSiteParsedList, self.atomTypeElements, self.structureResidueD = None, None, None
        self.remarksD, self.remarksAllD = None, None
        self.alignTarget, self.alignTemplate, self.alignOffsets = None, None, None
        self.templateEntityId, self.targetEntityId, self.templateDataId, self.targetDataId, self.alignmentDataId, self.coordDataId = None, None, None, None, None, None

    def convertPdbToCif(self, pdbFileIn, cifFileOut, **kwargs):
        self.organism = kwargs.get("organism", self.organism)
        self.uniProtId = kwargs.get("uniProtId", self.uniProtId)
        readOk = self.readPdb(pdbFileIn=pdbFileIn, organism=self.organism, uniProtId=self.uniProtId)
        writeOk = self.writeCif(cifFileOut=cifFileOut, organism=self.organism, uniProtId=self.uniProtId)
        return bool(readOk and writeOk)

    def readPdb(self, pdbFileIn, **kwargs):
        self.organism = kwargs.get("organism", self.organism)
        self.uniProtId = kwargs.get("uniProtId", self.uniProtId)
        self.pdbRead = False
        self.title, self.expdta, self.author, self.revDate = None, None, None, None
        self.jrnl, self.remarks, self.atoms = [], [], []
        self.align = False
        self.todaysDate = datetime.datetime.today().strftime("%Y-%m-%d")
        self.modelFileName = self.__fU.getFileName(pdbFileIn)
        self.modelNameBase = self.modelFileName.split('.pdb')[0]
        #
        try:
            with open(pdbFileIn, "r", encoding="utf-8") as fh:
                for line in fh:
                    # Handle standard SWISS-MODEL headers
                    if line.startswith('TITLE     '):
                        self.title = line[10:].strip()
                    elif line.startswith('EXPDTA    '):
                        self.expdta = line[10:].strip()
                    elif line.startswith('AUTHOR    '):
                        self.author = line[10:].strip()
                    elif line.startswith("REVDAT"):
                        revDate = line[10:].strip().split()[0]
                        self.revDate = datetime.datetime.strptime(revDate, "%d-%b-%y").strftime("%Y-%m-%d")
                        # self.pdbRevInfo = line[10:].strip()  # Then process this info separately to 'database_PDB_rev' category items
                    elif line.startswith('JRNL     '):
                        self.jrnl.append(line.strip())
                        # Then process these separately to get References, Disclaimers, etc.
                        # For now, just print out pre-organized citation information (assuming citation details won't change for a while)
                    elif line.startswith('REMARK'):
                        self.remarks.append(line.strip())
                    elif line.startswith('ATOM') or line.startswith('HETATM'):
                        self.atoms.append(line)
            # NOT ALL SWISS-MODEL models are single chain
            self.atomSiteParsedList, self.atomTypeElements = self.__parseAtomSites()
            self.structureResidueD = self.__getStructureResidueInfo()
            self.remarksD, self.remarksAllD = self.__parseRemarks()
            self.alignTarget, self.alignTemplate, self.alignOffsets = self.__parseAlign()
            assert len(self.alignTemplate) > 0
            self.align = True
            self.templateEntityId, self.targetEntityId = 1, 2
            self.templateDataId, self.targetDataId = 1, 2
            self.alignmentDataId, self.coordDataId = 3, 4
            self.pdbRead = True
        #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return self.pdbRead

    def __parseRemarks(self):
        remarksD, remarksAllD, remarkIdx = {}, {}, 0
        collectInfoFlag = False
        for ln in self.remarks:
            idx = ln.split()[1]
            if idx != remarkIdx:
                remarkIdx = idx
                remarksAllD[remarkIdx] = ""
            else:
                remarksAllD[remarkIdx] += ln[11:] + "\n"
            if len(ln) <= 11:
                continue
            if ln[10:].strip() == "MODEL INFORMATION":
                collectInfoFlag = True
                infoKey = "MODEL INFORMATION"
                remarksD[infoKey] = {}
                continue
            if ln[10:20].strip() == "TEMPLATE":
                collectInfoFlag = True
                infoKey = "TEMPLATE"
                # infoKey = ln[10:].strip()
                remarksD[infoKey] = {}
                continue
            if ln[10:15].strip() == "ALN" and ln[18:21] == "TRG":
                alnTrgKey = ln[10:21].strip()
                if alnTrgKey not in remarksD:
                    collectInfoFlag = False
                    remarksD[alnTrgKey] = ""
                remarksD[alnTrgKey] += ln[21:].strip()
            if ln[10:15].strip() == "ALN" and ln[18:21] == "TPL":
                alnTplKey = ln[10:21].strip()
                if alnTplKey not in remarksD:
                    collectInfoFlag = False
                    remarksD[alnTplKey] = ""
                remarksD[alnTplKey] += ln[21:].strip()
            if ln[10:15].strip() == "ALN" and ln[18:21] == "OFF":
                alnOffKey = ln[10:21].strip()
                remarksD[alnOffKey] = ln[21:].strip()
            if collectInfoFlag:
                key, val = ln[10:20].strip(), ln[20:].strip()
                remarksD[infoKey][key] = val
        return remarksD, remarksAllD

    def __parseAlign(self):
        # Search through ALN remarks and deterine whether Target or Template, and the Chain IDs A, B, etc...
        for k in self.remarksD:
            if k.startswith("ALN"):
                alignType = k[-3:]
                alignChain = k[4]
                if alignType == "TRG":
                    targetAlignD = {"chainId": alignChain, "oneLetterWithGaps": self.remarksD[k], "oneLetter": self.remarksD[k].replace('-', '')}
                if alignType == "TPL":
                    templateAlignD = {"chainId": alignChain, "oneLetterWithGaps": self.remarksD[k], "oneLetter": self.remarksD[k].replace('-', '')}
                if alignType == "OFF":
                    offAlignD = {"chainId": alignChain, "offset": self.remarksD[k]}
        return targetAlignD, templateAlignD, offAlignD

    def __parseAtomSites(self):
        atomSiteParsedList = []
        chainIds = [letter for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"]
        atomTypeElements = set()
        entityId = self.targetEntityId
        resNum = 1
        chainIdx = 0
        ordinal = 1
        authResNum = None
        authChainId = None
        for atom in self.atoms:
            # Detect new residue if PDB resnum changed
            thisResNum = atom[22:26].strip()
            thisChainId = atom[21].strip()
            if authResNum is not None and thisResNum != authResNum:
                resNum += 1
            # Detect new chain if chainID changed (and reset resNum back to 1)
            if authChainId is not None and thisChainId != authChainId:
                chainIdx += 1
                resNum = 1
            authResNum = thisResNum
            authChainId = thisChainId
            insCode = atom[26:27].strip() or '?'
            element = atom[76:78].strip()
            atomTypeElements.add(element)
            aD = {
                "group_PDB": atom[:6].strip(),
                "id": ordinal,
                "type_symbol": element,
                "label_atom_id": atom[12:16],
                "label_alt_id": ".",
                "label_comp_id": atom[17:20],
                "label_asym_id": chainIds[chainIdx],
                "label_entity_id": entityId,
                "label_seq_id": resNum,
                "pdbx_PDB_ins_code": insCode,
                "Cartn_x": atom[30:38],
                "Cartn_y": atom[38:46],
                "Cartn_z": atom[46:54],
                "occupancy": atom[54:60],
                "B_iso_or_equiv": atom[60:66],
                "auth_seq_id": authResNum,
                "auth_comp_id": atom[17:20],
                "auth_asym_id": authChainId,
            }
            atomSiteParsedList.append(aD)
            ordinal += 1
        return atomSiteParsedList, atomTypeElements

    def __getStructureResidueInfo(self):
        """Get residue-level details for the current structure (using the coordinates)"""
        pdbSeqIds, pdbChainIds, authSeqIds, authChainIds, pdbRes3Letter, pdbRes1Letter, insCodes, bFactors = [], [], [], [], [], [], [], []
        resTracker = None
        # Loop over atoms to get residue-level details
        for aD in self.atomSiteParsedList:
            if aD["auth_seq_id"] != resTracker and aD["group_PDB"] == "ATOM":
                resTracker = aD["auth_seq_id"]
                pdbSeqIds.append(aD["label_seq_id"])
                pdbChainIds.append(aD["label_asym_id"])
                authSeqIds.append(aD["auth_seq_id"])
                authChainIds.append(aD["auth_asym_id"])
                pdbRes3Letter.append(aD["label_comp_id"])
                pdbRes1Letter.append(three_to_one[aD["label_comp_id"]])
                bFactors.append(aD["B_iso_or_equiv"])
                if aD["pdbx_PDB_ins_code"] == "?":
                    insCodes.append(".")
                else:
                    insCodes.append(aD["pdbx_PDB_ins_code"])
        structureResidueD = {"pdbSeqIds": pdbSeqIds, "pdbChainIds": pdbChainIds, "authSeqIds": authSeqIds, "authChainIds": authChainIds,
                             "pdbRes3Letter": pdbRes3Letter, "pdbRes1Letter": pdbRes1Letter, "insCodes": insCodes, "bFactors": bFactors}
        return structureResidueD

    def print(self, strP):
        print(strP, file=self.fh)

    def loop(self, category, keys):
        return CifLoop(self.fh, category, keys)

    def writeCif(self, cifFileOut, **kwargs):
        try:
            self.organism = kwargs.get("organism", self.organism)
            self.uniProtId = kwargs.get("uniProtId", self.uniProtId)
            if not self.pdbRead:
                return False
            if not self.organism:
                self.organism = "?"
            with open(cifFileOut, "w", encoding="utf-8") as fh:
                self.fh = fh
                self.__writeHeader()
                self.__writeExptl()
                self.__writeCitationDetails()
                self.__writeChemComp()
                self.__writeEntityDetails()
                self.__writeTemplateDetails()
                self.__writeTargetDetails()
                self.__writeAlignment()
                self.__writeAssembly()
                self.__writeData()
                self.__writeSoftware()
                self.__writePdbRemarkText()
                self.__writeModelList()  # Check if this is writing what it should, and below
                self.__writeScores()
                self.__writeTargetRefDbDetails()
                self.__writePdbxDatabaseStatus()
                self.__writePdbxAuditRevisionDetails()
                self.__writePdbxAuditRevisionHistory()
                self.__writeAsym()
                self.__writePdbxPolySeqScheme()
                self.__writeAtomSites()
                fh.close()
                return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            return False

    def __writeHeader(self):
        self.print(f"data_{self.modelNameBase}")
        self.print(f"#\n_entry.id {self.modelNameBase}")
        self.print(f"_struct.entry_id {self.modelNameBase}")
        if self.title:
            self.print(f"_struct.title '{self.title}'")

    def __writeExptl(self):
        if self.expdta.startswith('THEORETICAL MODEL'):
            self.print(f"#\n_exptl.entry_id {self.modelNameBase}")
            self.print("_exptl.method 'THEORETICAL MODEL'")
            self.print(f"_exptl.details '{self.expdta}'")

    def __writeCitationDetails(self):
        """Citation details originally generated by MAXIT using a SWISS-MODEL PDB file created on 13-JUL-2021."""
        self.fh.write("""#
_audit_author.name           "SWISS-MODEL SERVER (SEE REFERENCE IN JRNL Records)"
_audit_author.pdbx_ordinal   1
#
loop_
_citation.id
_citation.title
_citation.journal_abbrev
_citation.journal_volume
_citation.page_first
_citation.page_last
_citation.year
_citation.journal_id_ASTM
_citation.country
_citation.journal_id_ISSN
_citation.journal_id_CSD
_citation.book_publisher
_citation.pdbx_database_id_PubMed
_citation.pdbx_database_id_DOI
primary "Swiss-Model: Homology Modelling Of Protein Structures And Complexes"                                             \
Nucleic.Acids.Res.. 46 W296 ? 2018 NARHAD UK 0305-1048 0389 ? 29788355 10.1093/nar/gky427
1       "The Swiss-Model Repository - New Features And Functionality"                                                     \
Nucleic.Acids.Res.. 45 ?    ? 2017 NARHAD UK 0305-1048 0389 ? 27899672 10.1093/nar/gkw1132
2       "Automated Comparative Protein Structure Modeling With Swiss-Model And Swiss-Pdbviewer: A Historical Perspective" \
Electrophoresis     30 ?    ? 2009 ?      GE 0173-0835 ?    ? 19517507 10.1002/elps.200900140
3       "Promod3 - A Versatile Homology Modelling Toolbox"                                                                \
"Plos Comp. Biol."  17 ?    ? 2021 ?      ?  ?         ?    ? 33507980 https://doi.org/10.1371/journal.pcbi.1008667
4       "Qmeandisco - Distance Constraints Applied On Model Quality Estimation"                                           \
Bioinformatics      36 ?    ? 2020 ?      UK 1460-2059 ?    ? 31697312 https://doi.org/10.1093/bioinformatics/btz828
5       "Toward The Estimation Of The Absolute Quality Of Individual Protein Structure Models"                            \
Bioinformatics      27 ?    ? 2011 ?      UK 1367-4803 ?    ? 21134891 10.1093/bioinformatics/btq662
6       "Modeling Protein Quaternary Structure Of Homo- And Hetero-Oligomers Beyond Binary Interactions By Homology"      \
Sci.Rep.            7  ?    ? 2017 ?      UK 2045-2322 ?    ? 28874689 10.1038/s41598-017-09654-8
#
loop_
_citation_author.citation_id
_citation_author.name
_citation_author.identifier_ORCID
_citation_author.ordinal
primary "Waterhouse, A."  ? 1
primary "Bertoni, M."     ? 2
primary "Bienert, S."     ? 3
primary "Studer, G."      ? 4
primary "Tauriello, G."   ? 5
primary "Gumienny, R."    ? 6
primary "Heer, F.T."      ? 7
primary "De Beer, T.A.P." ? 8
primary "Rempfer, C."     ? 9
primary "Bordoli, L."     ? 10
primary "Lepore, R."      ? 11
primary "Schwede, T."     ? 12
1       "Bienert, S."     ? 13
1       "Waterhouse, A."  ? 14
1       "De Beer, T.A.P." ? 15
1       "Tauriello, G."   ? 16
1       "Studer, G."      ? 17
1       "Bordoli, L."     ? 18
1       "Schwede, T."     ? 19
2       "Guex, N."        ? 20
2       "Peitsch, M.C."   ? 21
2       "Schwede, T."     ? 22
3       "Studer, G."      ? 23
3       "Tauriello, G."   ? 24
3       "Bienert, S."     ? 25
3       "Biasini, M."     ? 26
3       "Johner, N."      ? 27
3       "Schwede, T."     ? 28
4       "Studer, G."      ? 29
4       "Rempfer, C."     ? 30
4       "Waterhouse, A."  ? 31
4       "Gumienny, R."    ? 32
4       "Haas, J."        ? 33
4       "Schwede, T."     ? 34
5       "Benkert, P."     ? 35
5       "Biasini, M."     ? 36
5       "Schwede, T."     ? 37
6       "Bertoni, M."     ? 38
6       "Kiefer, F."      ? 39
6       "Biasini, M."     ? 40
6       "Bordoli, L."     ? 41
6       "Schwede, T."     ? 42
#""")

    def __writeChemComp(self):
        # Assume all 20 standard amino acids are in the model
        with self.loop(
                "chem_comp",
                ["id", "type", "name", "formula", "formula_weight"]) as lp:
            lp.write("""ALA 'L-peptide linking' ALANINE 'C3 H7 N O2' 89.094
ARG 'L-peptide linking' ARGININE 'C6 H15 N4 O2 1' 175.212
ASN 'L-peptide linking' ASPARAGINE 'C4 H8 N2 O3' 132.119
ASP 'L-peptide linking' 'ASPARTIC ACID' 'C4 H7 N O4' 133.103
CYS 'L-peptide linking' CYSTEINE 'C3 H7 N O2 S' 121.154
GLN 'L-peptide linking' GLUTAMINE 'C5 H10 N2 O3' 146.146
GLU 'L-peptide linking' 'GLUTAMIC ACID' 'C5 H9 N O4' 147.130
GLY 'peptide linking' GLYCINE 'C2 H5 N O2' 75.067
HIS 'L-peptide linking' HISTIDINE 'C6 H10 N3 O2 1' 156.165
ILE 'L-peptide linking' ISOLEUCINE 'C6 H13 N O2' 131.175
LEU 'L-peptide linking' LEUCINE 'C6 H13 N O2' 131.175
LYS 'L-peptide linking' LYSINE 'C6 H15 N2 O2 1' 147.198
MET 'L-peptide linking' METHIONINE 'C5 H11 N O2 S' 149.208
PHE 'L-peptide linking' PHENYLALANINE 'C9 H11 N O2' 165.192
PRO 'L-peptide linking' PROLINE 'C5 H9 N O2' 115.132
SER 'L-peptide linking' SERINE 'C3 H7 N O3' 105.093
THR 'L-peptide linking' THREONINE 'C4 H9 N O3' 119.120
TRP 'L-peptide linking' TRYPTOPHAN 'C11 H12 N2 O2' 204.229
TYR 'L-peptide linking' TYROSINE 'C9 H11 N O3' 181.191
VAL 'L-peptide linking' VALINE 'C5 H11 N O2' 117.148""")

    def __writeEntityDetails(self):
        with self.loop(
                "entity",
                ["id", "type", "src_method", "pdbx_description"]) as lp:
            if self.align:
                lp.write(f"{self.templateEntityId} polymer man template")
            lp.write(f"{self.targetEntityId} polymer man target")

        sD = self.structureResidueD
        targetPrimary = "".join(sD["pdbRes1Letter"][i] for i in range(len(sD["pdbRes1Letter"])) if sD["pdbChainIds"][i] == "A")

        with self.loop(
                "entity_poly",
                ["entity_id", "type", "nstd_linkage", "pdbx_seq_one_letter_code", "pdbx_seq_one_letter_code_can"]) as lp:
            if self.align:
                if self.alignTarget["oneLetter"].find(targetPrimary) < 0:
                    raise ValueError("Model sequence does not match target sequence in alignment:", targetPrimary, self.alignTarget["oneLetter"])
                pSeq = self.alignTemplate["oneLetter"]
                lp.write(f"{self.templateEntityId} polypeptide(L) no {pSeq} {pSeq}")
                lp.write(f"{self.targetEntityId} polypeptide(L) no {targetPrimary} {targetPrimary}")

        with self.loop(
                "entity_poly_seq",
                ["entity_id", "num", "mon_id", "hetero"]) as lp:
            if self.align:
                pSeq = self.alignTemplate["oneLetter"]
                for i, val in enumerate(pSeq):
                    lp.write(f"{self.templateEntityId} {i+1} {one_to_three[val]} n")
            for i, val in enumerate(self.structureResidueD["pdbRes3Letter"]):
                if self.structureResidueD["pdbChainIds"][i] == "A":  # Only write out chain A
                    lp.write(f'{self.targetEntityId} {self.structureResidueD["pdbSeqIds"][i]} {val} n')

    def __writePdbxPolySeqScheme(self):
        sD = self.structureResidueD
        with self.loop(
                "pdbx_poly_seq_scheme",
                ["asym_id", "entity_id", "seq_id", "mon_id", "ndb_seq_num", "pdb_seq_num", "auth_seq_num", "pdb_mon_id", "auth_mon_id",
                 "pdb_strand_id", "pdb_ins_code", "hetero"]) as lp:
            for i, val in enumerate(sD["pdbRes3Letter"]):
                lp.write("%s %d %d %s %d %s %s %s %s %s %s n" %
                         (sD["pdbChainIds"][i], self.targetEntityId, sD["pdbSeqIds"][i], val, sD["pdbSeqIds"][i],
                          sD["authSeqIds"][i], sD["authSeqIds"][i], val, val, sD["authChainIds"][i], sD["insCodes"][i]))

    def __writeTemplateDetails(self):
        if not self.align:
            return
        # Define the identity transformation (id=1)
        self.print("#\n_ma_template_trans_matrix.id 1")
        for i in range(1, 4):
            for j in range(1, 4):
                self.print("_ma_template_trans_matrix.rot_matrix[%d][%d] %s" % (j, i, "1.0" if i == j else "0.0"))
        for i in range(1, 4):
            self.print("_ma_template_trans_matrix.tr_vector[%d] 0.0" % i)

        with self.loop(
                "ma_template_details",
                ["ordinal_id", "template_id", "template_origin", "template_entity_type", "template_trans_matrix_id",
                 "template_data_id", "target_asym_id", "template_label_asym_id", "template_label_entity_id",
                 "template_model_num"]) as lp:
            # template structure is data_id=1
            # trans_matrix_id=1 is the identity transformation
            # model_num=1
            lp.write('1 1 "reference database" polymer 1 %d %s %s 1 1' % (self.templateDataId, self.alignTarget["chainId"], self.alignTemplate["chainId"]))

        with self.loop(
                "ma_template_poly", ["template_id", "seq_one_letter_code", "seq_one_letter_code_can"]) as lp:
            pSeq = self.alignTemplate["oneLetter"]
            lp.write("1 %s %s" % (pSeq, pSeq))

        with self.loop(
                "ma_template_poly_segment",
                ["id", "template_id"]) as lp:
            lp.write("1 1")

        with self.loop(
                "ma_template_ref_db_details",
                ["template_id", "db_name", "db_accession_code"]) as lp:
            lp.write("1 PDB %s" % self.remarksD["TEMPLATE"]["PDBID"])

    def __writeTargetDetails(self):
        with self.loop(
                "ma_target_entity", ["entity_id", "data_id", "origin"]) as lp:
            lp.write("%d %d ." % (self.targetEntityId, self.targetDataId))

        with self.loop(
                "ma_target_entity_instance",
                ["asym_id", "entity_id", "details"]) as lp:
            lp.write("%s %d ." % (self.alignTarget["chainId"], self.targetEntityId))

        if self.align:
            # Cannot write a template segment ID without an alignment
            with self.loop(
                    "ma_target_template_poly_mapping",
                    ["id", "template_segment_id", "target_asym_id", "target_seq_id_begin", "target_seq_id_end"]) as lp:
                lp.write("1 1 %s 1 %d" % (self.alignTarget["chainId"], len(self.structureResidueD["pdbRes1Letter"])))

    def __writeAlignment(self):
        if not self.align:
            return
        with self.loop(
                "ma_alignment_info",
                ["alignment_id", "data_id", "software_id", "alignment_length", "alignment_type", "alignment_mode"]) as lp:
            lp.write('1 %d 1 %d "target-template pairwise alignment" global'
                     % (self.alignmentDataId, len(self.alignTemplate["oneLetterWithGaps"])))

        with self.loop(
                "ma_alignment_details",
                # Add in "score_value" after determining which SwissModel header field represents the 'e-value'
                ["ordinal_id", "alignment_id", "template_segment_id", "target_asym_id", "score_type", "percent_sequence_identity"]) as lp:
            lp.write("1 1 1 %s '%s e-value' %s" % (self.alignTarget["chainId"], self.remarksD["TEMPLATE"]["FOUND"], self.remarksD["TEMPLATE"]["SID"]))

        with self.loop(
                "ma_alignment",
                ["ordinal_id", "alignment_id", "target_template_flag", "sequence"]) as lp:
            # Template (flag=2)
            lp.write("1 1 2 %s" % self.alignTemplate["oneLetterWithGaps"])
            # Target (flag=1)
            lp.write("2 1 1 %s" % self.alignTarget["oneLetterWithGaps"])

    def __writeAssembly(self):
        # This may need to be looped over for each chain, as done in 'struct_asym' below
        with self.loop(
                'ma_struct_assembly',
                ['ordinal_id', 'assembly_id', 'entity_id', 'asym_id', 'seq_id_begin', 'seq_id_end']) as lp:
            # Simple assembly of a single chain
            lp.write("1 1 %d %s 1 %d" % (self.targetEntityId, self.alignTarget["chainId"], len(self.structureResidueD["pdbRes1Letter"])))

    def __writeData(self):
        with self.loop("ma_data", ["id", "name", "content_type"]) as lp:
            lp.write("%d 'Template Structure' 'template structure'" % self.templateDataId)
            lp.write("%d 'Target Sequence' target" % self.targetDataId)
            lp.write("%d 'Target Template Alignment' 'target-template alignment'" % self.alignmentDataId)
            lp.write("%d 'Target Structure' 'model coordinates'" % self.coordDataId)

        # Put each data item in its own group
        with self.loop("ma_data_group", ["ordinal_id", "group_id", "data_id"]) as lp:
            for i in range(1, 5):
                lp.write(f"{i} {i} {i}")

    def __writeSoftware(self):
        with self.loop(
                "software",
                ["pdbx_ordinal", "name", "classification", "version", "type", "location", "citation_id"]) as lp:
            lp.write("1 SWISS-MODEL 'comparative modeling' %s program https://swissmodel.expasy.org/ primary"
                     % self.remarksD["TEMPLATE"]["SMTLV"])  # Version of SWISS-MODEL Template Library
            lp.write("2 ProMod3 'comparative modeling' %s program https://openstructure.org/promod3/3.2/ 3"
                     % self.remarksD["MODEL INFORMATION"]["VERSN"])

        # Put each piece of software in its own group
        with self.loop(
                "ma_software_group",
                ["ordinal_id", "group_id", "software_id"]) as lp:
            for i in range(1, 3):
                lp.write(f"{i} {i} {i}")

    def __writeScores(self):
        with self.loop(
                'ma_qa_metric',
                ['id', 'name', 'description', 'type', 'mode', 'type_other_details', 'software_group_id']) as lp:
            lp.write("1 GMQE 'Global Model Quality Estimate' other global . 1")  # SWISS-MODEL is software_id=1
            lp.write("2 QMEANDisCo 'Average Model Confidence' other global 'QMEAN version %s' 1" % self.remarksD["MODEL INFORMATION"]["QMNV"])  # This is "QMNDG" in remarks
            lp.write("3 QMEANDisCo 'Per-residue Quality Score' other local 'QMEAN version %s' 1" % self.remarksD["MODEL INFORMATION"]["QMNV"])  # This is beta factor field
        with self.loop(
                'ma_qa_metric_global',
                ['ordinal_id', 'model_id', 'metric_id', 'value']) as lp:
            lp.write("1 1 1 %s" % self.remarksD["MODEL INFORMATION"]["GMQE"])
            lp.write("2 1 2 %s" % self.remarksD["MODEL INFORMATION"]["QMNDG"])
        with self.loop(
                'ma_qa_metric_local',
                ['ordinal_id', 'model_id', 'label_asym_id', 'label_comp_id', 'label_seq_id', 'metric_id', 'metric_value']) as lp:
            sD = self.structureResidueD
            for resIdx in range(len(sD["pdbRes1Letter"])):
                # Per-residue quality score is metric_id 3
                lp.write("1 1 %s %s %d 3 %s" % (sD["pdbChainIds"][resIdx], sD["pdbRes3Letter"][resIdx], sD["pdbSeqIds"][resIdx], sD["bFactors"][resIdx]))

    def __writeModelList(self):
        with self.loop(
                'ma_model_list',
                ['ordinal_id', 'model_id', 'model_group_id', 'model_name', 'model_group_name', 'assembly_id', 'data_id', 'model_type']) as lp:
            lp.write(f"1 1 1 'Selected model' . 1 {self.coordDataId} 'Homology model'")

    def __writeTargetRefDbDetails(self):
        if not self.uniProtId:
            return
        sD = self.structureResidueD
        self.print("#\n_ma_target_ref_db_details.db_name UNP")
        self.print(f"_ma_target_ref_db_details.accession {self.uniProtId}")
        self.print(f"_ma_target_ref_db_details.code {self.uniProtId}")  # To be consitent with AF this should be the "ID" value mapped from "ACC+ID" UniProt ID (or "LOCUS" on NCBI)
        if self.organism:
            self.print(f"_ma_target_ref_db_details.organism_scientific {self.organism}")
        self.print(f'_ma_target_ref_db_details.seq_db_align_begin {sD["authSeqIds"][0]}')
        self.print(f'_ma_target_ref_db_details.seq_db_align_end {sD["authSeqIds"][-1]}')
        self.print(f"_ma_target_ref_db_details.target_entity_id {self.targetEntityId}")

    def __writePdbxAuditRevisionDetails(self):
        self.print('#\n_pdbx_audit_revision_details.data_content_type "Structure model"')
        self.print('_pdbx_audit_revision_details.description ?')
        self.print('_pdbx_audit_revision_details.ordinal 1')
        self.print('_pdbx_audit_revision_details.provider repository')
        self.print('_pdbx_audit_revision_details.revision_ordinal 1')
        self.print('_pdbx_audit_revision_details.type "Initial release"')

    def __writePdbxAuditRevisionHistory(self):
        self.print('#\n_pdbx_audit_revision_history.data_content_type "Structure model"')
        self.print('_pdbx_audit_revision_history.major_revision 1')
        self.print('_pdbx_audit_revision_history.minor_revision 0')
        self.print('_pdbx_audit_revision_history.ordinal 1')
        self.print(f'_pdbx_audit_revision_history.revision_date {self.revDate}')

    def __writePdbxDatabaseStatus(self):
        self.print(f'#\n_pdbx_database_status.entry_id {self.targetEntityId}')
        self.print(f'_pdbx_database_status.recvd_initial_deposition_date {self.todaysDate}')
        self.print('_pdbx_database_status.status_code REL')

    def __writeAsym(self):
        # Revist '__writeAssembly' too
        sD = self.structureResidueD
        uniqPdbChainIds, uniqAuthChainIds = [], []
        for chain in sD["pdbChainIds"]:
            if chain not in uniqPdbChainIds:
                uniqPdbChainIds.append(chain)
        for chain in sD["authChainIds"]:
            if chain not in uniqAuthChainIds:
                uniqAuthChainIds.append(chain)
        assert len(uniqPdbChainIds) == len(uniqAuthChainIds)
        with self.loop('struct_asym', ['id', 'entity_id', 'pdbx_PDB_id', 'pdbx_alt_id', 'details']) as lp:
            for i in range(len(uniqPdbChainIds)):
                lp.write("%s %s %s %d ?" % (uniqPdbChainIds[i], uniqAuthChainIds[i], uniqAuthChainIds[i], self.targetEntityId))

    def __writeAtomSites(self):
        with self.loop(
            'atom_site',
            ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "label_asym_id",
             "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
             "B_iso_or_equiv", "auth_seq_id", "auth_asym_id", "auth_comp_id"]) as lp:
            for aD in self.atomSiteParsedList:
                lp.write("%s %d %s %s %s %s %s %d %d %s %s %s %s %s %s %s %s %s"
                         % (aD["group_PDB"], aD["id"], aD["type_symbol"], aD["label_atom_id"], aD["label_alt_id"], aD["label_comp_id"], aD["label_asym_id"],
                            aD["label_entity_id"], aD["label_seq_id"], aD["pdbx_PDB_ins_code"], aD["Cartn_x"], aD["Cartn_y"], aD["Cartn_z"], aD["occupancy"],
                            aD["B_iso_or_equiv"], aD["auth_seq_id"], aD["auth_asym_id"], aD["auth_comp_id"]))

        with self.loop('atom_type', ['symbol']) as lp:
            for element in sorted(self.atomTypeElements):
                lp.write(element)

    def __writePdbRemarkText(self):
        with self.loop("database_PDB_remark", ["id", "text"]) as lp:
            for rem in self.remarksAllD:
                lp.write("%s\n;%s;" % (rem, self.remarksAllD[rem]))
