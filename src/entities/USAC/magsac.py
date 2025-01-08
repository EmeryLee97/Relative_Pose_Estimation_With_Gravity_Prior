# requirments: OpenCV >= 4.0

import cv2
import numpy as np

def run_magsacpp(pts_1: np.ndarray, pts_2: np.ndarray, intrinsics=np.eye(3)):
    """
    Args:
        pts_1, pts_2: (n, 2) point sets
        intrinsics: (3, 3) camera intrinsic matrix

    Returns:
        relative rotation and translation
    """

    E, _ = cv2.findEssentialMat(
        pts_1, 
        pts_2, 
        intrinsics,
        # method=cv2.USAC_MAGSAC,
        method=cv2.USAC_ACCURATE,
        prob=0.99, 
        threshold=1.0
    )
    _, rot, trans, _ = cv2.recoverPose(E, pts_1, pts_2, intrinsics)

    return rot, trans

def run_gcransac(pts_1: np.ndarray, pts_2: np.ndarray, intrinsics=np.eye(3)):
    """
    Args:
        pts_1, pts_2: (n, 2) point sets
        intrinsics: (3, 3) camera intrinsic matrix

    Returns:
        relative rotation and translation
    """

    E, _ = cv2.findEssentialMat(
        pts_1, 
        pts_2, 
        intrinsics,
        method=cv2.USAC_ACCURATE,
        prob=0.99, 
        threshold=1.0
    )
    _, rot, trans, _ = cv2.recoverPose(E, pts_1, pts_2, intrinsics)

    return rot, trans