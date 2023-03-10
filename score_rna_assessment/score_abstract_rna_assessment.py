"""
Class that get and convert the pdb files to structures encoded by RNA-tools.
"""
from typing import Optional, Tuple

from RNA_normalizer.structures.pdb_struct import PDBStruct


class ScoreAbstractRnaAssessment:
    def __init__(self, *args, **kwargs):
        super(ScoreAbstractRnaAssessment).__init__(*args, **kwargs)

    @staticmethod
    def convert_pdb_to_structure(
        pred_path: str,
        native_path: str,
        native_index: Optional[str] = None,
        prediction_index: Optional[str] = None,
    ) -> Tuple[PDBStruct, PDBStruct]:
        """
        Convert the .pdb files to structures readable by RNA-tools.
        :param pred_path: the path to the .pdb file of a prediction.
        :param native_path: the path to the .pdb file of the native structure.
        :param native_index: file that describes the delimitation of the RNA for the native file
        :param prediction_index: file that describes the delimitation of the RNA
                    for the prediction file
        :return: two instances of PDBStruct for the native and prediction structures
        """
        native_struc, pred_struc = PDBStruct(), PDBStruct()
        native_struc.load(native_path, native_index)
        pred_struc.load(pred_path, prediction_index)
        return native_struc, pred_struc

