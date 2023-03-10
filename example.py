import os
from typing import Dict

from score_rna_assessment.score_di import ScoreDI
from score_rna_assessment.score_inf import ScoreINF
from score_rna_assessment.score_p_value import ScorePValue
from score_rna_assessment.score_rmsd import ScoreRMSD


def compute_all_score(pred_path: str, native_path: str) -> Dict:
    """
    Compute RMSD, P-value, INF scores and DI
    Args:
        :param pred_path: path to a .pdb file of a predicted structure
        :param native_path: path to a .pdb native structure
    Return:
        dictionary with the different scores for the given structures
    """
    rmsd = ScoreRMSD.compute_rmsd(pred_path, native_path)
    p_value = ScorePValue.compute_p_value(pred_path, native_path)
    inf_all, inf_wc, inf_nwc, inf_stack = (
        ScoreINF.compute_inf_all(pred_path, native_path),
        ScoreINF.compute_inf_wc(pred_path, native_path),
        ScoreINF.compute_inf_nwc(pred_path, native_path),
        ScoreINF.compute_inf_stack(pred_path, native_path),
    )
    di = ScoreDI.compute_di(pred_path, native_path)
    output = {
        "RMSD": rmsd,
        "P-VALUE": p_value,
        "INF_ALL": inf_all,
        "INF_WC": inf_wc,
        "INF_NWC": inf_nwc,
        "INF_STACK": inf_stack,
        "DI": di,
    }
    return output


if __name__ == "__main__":
    pred_path = os.path.join("example", "3dRNA-1Z43-10.pdb")
    native_path = os.path.join("example", "1Z43.pdb")
    scores = compute_all_score(pred_path, native_path)
    for name, score in scores.items():
        print(f"{name} : {score}")
