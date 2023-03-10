"""
Class that implements the Deformation Index score.
It uses the RNA_Assessment repo.
"""
from typing import  Optional

from RNA_normalizer.structures.pdb_struct import PDBStruct
from score_rna_assessment.score_abstract_rna_assessment import ScoreAbstractRnaAssessment
from score_rna_assessment.score_inf import ScoreINF
from score_rna_assessment.score_rmsd import ScoreRMSD


class ScoreDI(ScoreAbstractRnaAssessment):
    def __init__(self, *args, **kwargs):
        super(ScoreDI).__init__(*args, **kwargs)

    @staticmethod
    def compute_di_from_structures(native_struc: PDBStruct, pred_struc: PDBStruct) -> float:
        """
        Static method to compute the Deformation Index score
                from the native and predicted structures.
        This is defined as the RMSD / INF_all
        :param native_struc: native structure in a PDBStruc instance
        :param pred_struc: predicted structure in a PDBStruc instance
        :return: the DI score from these structures
        """
        rmsd = ScoreRMSD.compute_rmsd_from_structures(native_struc, pred_struc)
        inf_all = ScoreINF.compute_inf_all_from_structures(native_struc, pred_struc)
        return rmsd / inf_all

    @staticmethod
    def compute_di(
        pred_path: str,
        native_path: str,
        native_index: Optional[str] = None,
        prediction_index: Optional[str] = None,
    ) -> float:
        """
        Static method to compute the Deformation Index score
                from the native and predicted structures.
        This is defined as the RMSD / INF_all
        :param pred_path: the path to the .pdb file of a prediction.
        :param native_path: the path to the .pdb file of the native structure.
        :param native_index: file that describes the delimitation of the RNA for the native file
        :param prediction_index: file that describes the delimitation of the RNA
                    for the prediction file
        :return: the DI score from these structures
        """
        native_struc, pred_struc = ScoreAbstractRnaAssessment.convert_pdb_to_structure(
            pred_path, native_path, native_index, prediction_index
        )
        di = ScoreDI.compute_di_from_structures(native_struc, pred_struc)
        return di

