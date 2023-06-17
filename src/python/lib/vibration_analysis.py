""" Vibration analysis module.
"""
import numpy as np
from ase import units

speed_of_light = 137.035999084
au_to_kg = 9.1093837015e-31
amu_to_kg = 1.66053904e-27
amu_to_au = amu_to_kg / au_to_kg


# helper functions
def compute_moi(masses, pos_com_shifted):
    """Compute MOI (moment of inertia) tensor."""
    # Compute elements need to calculate MOI tensor
    mass_xyz_com_sq_sum = np.sum(masses[:, None] * pos_com_shifted**2, axis=0)
    mass_xy = np.sum(masses * pos_com_shifted[:, 0] * pos_com_shifted[:, 1], axis=0)
    mass_yz = np.sum(masses * pos_com_shifted[:, 1] * pos_com_shifted[:, 2], axis=0)
    mass_xz = np.sum(masses * pos_com_shifted[:, 0] * pos_com_shifted[:, 2], axis=0)

    # MOI tensor
    moi = np.array(
        [
            [
                mass_xyz_com_sq_sum[1] + mass_xyz_com_sq_sum[2],
                -1 * mass_xy,
                -1 * mass_xz,
            ],
            [
                -1 * mass_xy,
                mass_xyz_com_sq_sum[0] + mass_xyz_com_sq_sum[2],
                -1 * mass_yz,
            ],
            [
                -1 * mass_xz,
                -1 * mass_yz,
                mass_xyz_com_sq_sum[0] + mass_xyz_com_sq_sum[1],
            ],
        ]
    )
    return moi


def trans_rot_vec(masses, pos_com_shifted, moi_eigvec):
    # Mass-weighted translational vectors
    zero_vec = np.zeros([len(masses)])
    sqrtmassvec = np.sqrt(masses)
    expsqrtmassvec = np.repeat(sqrtmassvec, 3)

    d1 = np.transpose(np.stack((sqrtmassvec, zero_vec, zero_vec))).reshape(-1)
    d2 = np.transpose(np.stack((zero_vec, sqrtmassvec, zero_vec))).reshape(-1)
    d3 = np.transpose(np.stack((zero_vec, zero_vec, sqrtmassvec))).reshape(-1)

    # Mass-weighted rotational vectors
    P = np.matmul(pos_com_shifted, moi_eigvec)

    d4 = (
        np.repeat(P[:, 1], 3).reshape(-1)
        * np.tile(moi_eigvec[:, 2], len(masses)).reshape(-1)
        - np.repeat(P[:, 2], 3).reshape(-1)
        * np.tile(moi_eigvec[:, 1], len(masses)).reshape(-1)
    ) * expsqrtmassvec

    d5 = (
        np.repeat(P[:, 2], 3).reshape(-1)
        * np.tile(moi_eigvec[:, 0], len(masses)).reshape(-1)
        - np.repeat(P[:, 0], 3).reshape(-1)
        * np.tile(moi_eigvec[:, 2], len(masses)).reshape(-1)
    ) * expsqrtmassvec

    d6 = (
        np.repeat(P[:, 0], 3).reshape(-1)
        * np.tile(moi_eigvec[:, 1], len(masses)).reshape(-1)
        - np.repeat(P[:, 1], 3).reshape(-1)
        * np.tile(moi_eigvec[:, 0], len(masses)).reshape(-1)
    ) * expsqrtmassvec

    d1_norm = d1 / np.linalg.norm(d1)
    d2_norm = d2 / np.linalg.norm(d2)
    d3_norm = d3 / np.linalg.norm(d3)
    d4_norm = d4 / np.linalg.norm(d4)
    d5_norm = d5 / np.linalg.norm(d5)
    d6_norm = d6 / np.linalg.norm(d6)

    d_norms = np.stack((d1_norm, d2_norm, d3_norm, d4_norm, d5_norm, d6_norm))

    return d_norms


class VibAnalysis:
    def __init__(self, atoms, hessian):
        self.atoms = atoms
        self.hessian = hessian
        self.vib_freqs_cm = None
        self.vib_freqs_au = None
        self.hessian_eigvals = None
        self.hessian_eigvecs = None

    def _compute_vib_props(self):
        R = self.atoms.get_positions() / units.Bohr
        com = self.atoms.get_center_of_mass() / units.Bohr
        R_com = R - com
        M = self.atoms.get_masses() * amu_to_au
        M_3 = np.repeat(M, 3)
        H = self.hessian
        H_MWC = H / np.sqrt(np.multiply.outer(M_3, M_3))
        moi = compute_moi(M, R_com)
        _, moi_eigvecs = np.linalg.eig(moi)
        d_norms = trans_rot_vec(M, R_com, moi_eigvecs)

        P = np.identity(3 * len(self.atoms))
        for dx_norm in d_norms:
            P -= np.outer(dx_norm, dx_norm)

        H_MWC_proj = P.T @ H_MWC @ P

        H_eigvals, H_eigvecs = np.linalg.eigh(H_MWC_proj)
        H_eigval_abs = np.abs(H_eigvals)

        # displacements
        disp_cart = (
            np.einsum("ik,kj,i->ij", P, H_eigvecs, 1 / np.sqrt(M_3)) * units.Bohr
        )
        disp_cart /= np.linalg.norm(disp_cart, axis=0)
        self.displacements = []
        for i in range(len(self.atoms) * 3):
            self.displacements.append(disp_cart[:, i].reshape(-1, 3))
        self.displacements = np.stack(self.displacements)

        H_eigvals_neg_idx = np.where(H_eigvals < 0)[0]
        vib_freqs = np.sqrt(H_eigval_abs / (4 * np.pi**2 * speed_of_light**2))
        vib_freqs[H_eigvals_neg_idx] *= -1
        trans_rot_idx = np.where(np.abs(vib_freqs) < (units.Bohr / units.m) / 100)[0]
        vib_freqs = np.delete(vib_freqs, trans_rot_idx)
        force_constant = np.delete(H_eigvals, trans_rot_idx)
        self.displacements = np.delete(self.displacements, trans_rot_idx, axis=0)
        self.vib_freqs_cm = vib_freqs / (units.Bohr / units.m) / 100
        self.force_constant = force_constant
        self.hessian_eigvals = H_eigvals
        self.hessian_eigvecs = H_eigvecs
