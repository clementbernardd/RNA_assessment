"""
Class that computes the INF score.
This is based on the RNA-tools implementation.
The INF score was invented to assess the central characteristics of RNA architecture.
The INF can either measure different interaction types
    (WC base-pairing, non-WC base pairing, base stacking) separately or combine
    all of the types (resulting in INFwc, INFnwc, INFstacking and INFall)
"""
from typing import Optional

from RNA_normalizer.structures.pdb_comparer import PDBComparer
from RNA_normalizer.structures.pdb_struct import PDBStruct
from score_rna_assessment.score_abstract_rna_assessment import ScoreAbstractRnaAssessment


class ScoreINF(ScoreAbstractRnaAssessment):
    def __init__(self, *args, **kwargs):
        super(ScoreINF).__init__(*args, **kwargs)

    @staticmethod
    def compute_inf_all_from_structures(native_struc: PDBStruct, pred_struc: PDBStruct) -> float:
        """
        Compute the INF score combining all the types.
        :param native_struc: native structure in a PDBStruc instance
        :param pred_struc: predicted structure in a PDBStruc instance
        :return: the INF score of the associated molecules
        """
        comparer = PDBComparer()
        inf_all = comparer.INF(src_struct=pred_struc, trg_struct=native_struc, type="ALL")
        return inf_all

    @staticmethod
    def compute_inf_all(
        pred_path: str,
        native_path: str,
        native_index: Optional[str] = None,
        prediction_index: Optional[str] = None,
    ) -> float:
        """
        Compute the INF score combining all the types.
        :param pred_path: the path to the .pdb file of a prediction.
        :param native_path: the path to the .pdb file of the native structure.
        :param native_index: file that describes the delimitation of the RNA for the native file
        :param prediction_index: file that describes the delimitation of the RNA
                    for the prediction file
        :return: the INF score of the associated molecules
        """
        native_struc, pred_struc = ScoreAbstractRnaAssessment.convert_pdb_to_structure(
            pred_path, native_path, native_index, prediction_index
        )
        inf_all = ScoreINF.compute_inf_all_from_structures(native_struc, pred_struc)
        return inf_all

    @staticmethod
    def compute_inf_wc_from_structures(native_struc: PDBStruct, pred_struc: PDBStruct) -> float:
        """
        Compute the INF score only for the WC base-pairing
        :param native_struc: native structure in a PDBStruc instance
        :param pred_struc: predicted structure in a PDBStruc instance
        :return: the INFwc score of the associated molecules
        """
        comparer = PDBComparer()
        inf_wc = comparer.INF(src_struct=pred_struc, trg_struct=native_struc, type="PAIR_2D")
        return inf_wc

    @staticmethod
    def compute_inf_wc(
        pred_path: str,
        native_path: str,
        native_index: Optional[str] = None,
        prediction_index: Optional[str] = None,
    ) -> float:
        """
        Compute the INF score only for the WC base-pairing
        :param pred_path: the path to the .pdb file of a prediction.
        :param native_path: the path to the .pdb file of the native structure.
        :param native_index: file that describes the delimitation of the RNA for the native file
        :param prediction_index: file that describes the delimitation of the RNA
                    for the prediction file
        :return: the INFwc score of the associated molecules
        """
        native_struc, pred_struc = ScoreAbstractRnaAssessment.convert_pdb_to_structure(
            pred_path, native_path, native_index, prediction_index
        )
        inf_wc = ScoreINF.compute_inf_wc_from_structures(native_struc, pred_struc)
        return inf_wc

    @staticmethod
    def compute_inf_nwc_from_structures(native_struc: PDBStruct, pred_struc: PDBStruct) -> float:
        """
        Compute the INF score only for the non WC base-pairing
        :param native_struc: native structure in a PDBStruc instance
        :param pred_struc: predicted structure in a PDBStruc instance
        :return: the INFwc score of the associated molecules
        """
        comparer = PDBComparer()
        inf_nwc = comparer.INF(src_struct=pred_struc, trg_struct=native_struc, type="PAIR_3D")
        return inf_nwc

    @staticmethod
    def compute_inf_nwc(
        pred_path: str,
        native_path: str,
        native_index: Optional[str] = None,
        prediction_index: Optional[str] = None,
    ) -> float:
        """
        Compute the INF score only for the non WC base-pairing
        :param pred_path: the path to the .pdb file of a prediction.
        :param native_path: the path to the .pdb file of the native structure.
        :param native_index: file that describes the delimitation of the RNA for the native file
        :param prediction_index: file that describes the delimitation of the RNA
                    for the prediction file
        :return: the INFwc score of the associated molecules
        """
        native_struc, pred_struc = ScoreAbstractRnaAssessment.convert_pdb_to_structure(
            pred_path, native_path, native_index, prediction_index
        )
        inf_nwc = ScoreINF.compute_inf_nwc_from_structures(native_struc, pred_struc)
        return inf_nwc

    @staticmethod
    def compute_inf_stack_from_structures(native_struc: PDBStruct, pred_struc: PDBStruct) -> float:
        """
        Compute the INF score only for the non WC base-pairing
        :param native_struc: native structure in a PDBStruc instance
        :param pred_struc: predicted structure in a PDBStruc instance
        :return: the INFwc score of the associated molecules
        """
        comparer = PDBComparer()
        inf_stack = comparer.INF(src_struct=pred_struc, trg_struct=native_struc, type="STACK")
        return inf_stack

    @staticmethod
    def compute_inf_stack(
        pred_path: str,
        native_path: str,
        native_index: Optional[str] = None,
        prediction_index: Optional[str] = None,
    ) -> float:
        """
        Compute the INF score only for the non WC base-pairing
        :param pred_path: the path to the .pdb file of a prediction.
        :param native_path: the path to the .pdb file of the native structure.
        :param native_index: file that describes the delimitation of the RNA for the native file
        :param prediction_index: file that describes the delimitation of the RNA
                    for the prediction file
        :return: the INFwc score of the associated molecules
        """
        native_struc, pred_struc = ScoreAbstractRnaAssessment.convert_pdb_to_structure(
            pred_path, native_path, native_index, prediction_index
        )
        inf_stack = ScoreINF.compute_inf_stack_from_structures(native_struc, pred_struc)
        return inf_stack

