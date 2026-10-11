"""금속배선 모델.

- 4-point probe 면저항 측정: 무한 박막, 유한 원형 웨이퍼(영상법), 유한차분 수치해
- 배선 RC 지연 (Plummer Ch. 11, Eq. 1)과 전류 밀도

단위는 함수마다 표기한다.
"""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

EPS0 = 8.8541878128e-12  # F/m

# 비저항 [uOhm cm] (Plummer Ch. 11)
RESISTIVITY_UOHM_CM = {"poly-Si": 500.0, "WSi2": 30.0, "W": 10.0, "Al": 3.0, "Cu": 2.2}


# ---------------------------------------------------------------------------
# 4-point probe
# ---------------------------------------------------------------------------
INFINITE_SHEET_CF = np.pi / np.log(2.0)   # 4.532


def _probe_positions(spacing, x0=0.0):
    return x0 + spacing * np.array([-1.5, -0.5, 0.5, 1.5]) + 0j


def four_point_probe_cf_disk(diameter, spacing, x0=0.0):
    """원형 웨이퍼 위 일직선 4-point probe의 보정 계수 CF = Rs / (V/I).

    탐침이 지름 방향으로 놓이고 탐침 중심이 웨이퍼 중심에서 x0만큼 떨어져 있다.
    가장자리는 전류가 빠져나가지 않는 경계이므로, 원 바깥의 반전점 a^2 / conj(z)에
    같은 부호의 영상 전류원을 두어 경계 조건을 맞춘다.
    """
    a = diameter / 2
    z1, z2, z3, z4 = _probe_positions(spacing, x0)
    img = lambda z: a**2 / np.conj(z) if z != 0 else np.inf

    def u(z, zs):                                           # 전류원 zs와 영상이 z에 만드는 퍼텐셜 (단위: I Rs / 2 pi)
        return -np.log(abs(z - zs)) - np.log(abs(z - img(zs)))

    bracket = (u(z2, z1) - u(z2, z4)) - (u(z3, z1) - u(z3, z4))
    return 2 * np.pi / bracket


def smits_center_cf(diameter, spacing):
    """Smits(1958)의 원형 시료 중앙 측정 보정 계수 (교재와 측정 표준에서 쓰는 닫힌 식)."""
    r2 = (diameter / spacing) ** 2
    return np.pi / (np.log(2.0) + np.log((r2 + 3) / (r2 - 3)))


def four_point_probe_fd(diameter, spacing, x0=0.0, nodes_per_spacing=10):
    """유한차분으로 원형 박막의 전위를 풀어 보정 계수를 구한다 (검증용).

    격자 간격 h에서 이웃 격자점 사이의 컨덕턴스는 1/Rs로 h와 무관하다(2D 박막).
    원 밖의 이웃은 전류가 흐르지 않는 경계로 보고 연결하지 않는다.
    """
    if nodes_per_spacing % 2:
        raise ValueError("탐침 위치(+/-0.5 s, +/-1.5 s)가 격자점에 놓이도록 nodes_per_spacing은 짝수여야 합니다.")
    h = spacing / nodes_per_spacing
    a = diameter / 2
    n = int(np.ceil(a / h)) + 2
    ix = np.arange(-n, n + 1)
    X, Y = np.meshgrid(ix * h, ix * h, indexing="ij")
    inside = X**2 + Y**2 <= a**2
    idx = -np.ones(X.shape, dtype=int)
    idx[inside] = np.arange(inside.sum())
    rows, cols, vals = [], [], []
    diag = np.zeros(inside.sum())
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        nb = np.roll(np.roll(idx, -di, axis=0), -dj, axis=1)
        ok = inside & (nb >= 0)
        i_, j_ = idx[ok], nb[ok]
        rows.append(i_)
        cols.append(j_)
        vals.append(-np.ones_like(i_, dtype=float))
        np.add.at(diag, i_, 1.0)
    rows.append(np.arange(len(diag)))
    cols.append(np.arange(len(diag)))
    vals.append(diag)
    G = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(len(diag), len(diag))).tocsr()

    def node(x):
        i = int(round(x / h)) + n
        return idx[i, n]

    pz = _probe_positions(spacing, x0).real
    k1, k2, k3, k4 = (node(p) for p in pz)
    b = np.zeros(len(diag))
    b[k1], b[k4] = 1.0, -1.0                                # I = 1, Rs = 1
    G = G.tolil()
    G[k4, :] = 0                                            # 기준 전위 고정 (특이 행렬 해소)
    G[k4, k4] = 1.0
    b[k4] = 0.0
    phi = spsolve(G.tocsr(), b)
    return 1.0 / (phi[k2] - phi[k3])


def resistivity_from_sheet(Rs_ohm_sq, thickness_nm):
    """비저항 [uOhm cm] = Rs x t."""
    return Rs_ohm_sq * thickness_nm * 1e-7 * 1e6


def sheet_from_resistivity(rho_uohm_cm, thickness_nm):
    """면저항 [Ohm/sq] = rho / t."""
    return rho_uohm_cm * 1e-6 / (thickness_nm * 1e-7)


# ---------------------------------------------------------------------------
# 배선 RC 지연과 전류 밀도
# ---------------------------------------------------------------------------
def rc_delay(rho_uohm_cm, length_mm, W_um, H_um, x_ox_um, L_s_um, K_ox, K_I=1.0):
    """전역 배선의 지연 시간 tau = 0.89 K_I K_ox eps0 rho L^2 (1/(H x_ox) + 1/(W L_s)) [s].

    배선 아래 기판과의 용량(두께 x_ox)과 옆 배선과의 용량(간격 L_s)을 함께 고려한다.
    K_I는 가장자리 전계(fringing) 보정 계수다.
    """
    rho = rho_uohm_cm * 1e-8                                # Ohm m
    L = np.asarray(length_mm, dtype=float) * 1e-3
    W, H, xo, Ls = (v * 1e-6 for v in (W_um, H_um, x_ox_um, L_s_um))
    return 0.89 * K_I * K_ox * EPS0 * rho * L**2 * (1.0 / (H * xo) + 1.0 / (W * Ls))


def line_resistance_capacitance(rho_uohm_cm, length_mm, W_um, H_um, x_ox_um, L_s_um, K_ox, K_I=1.0):
    """배선 저항 R [Ohm]과 용량 C [F]. 아래 판과 양옆 중 한쪽 배선과의 평행판 근사."""
    rho = rho_uohm_cm * 1e-8
    L = length_mm * 1e-3
    W, H, xo, Ls = (v * 1e-6 for v in (W_um, H_um, x_ox_um, L_s_um))
    R = rho * L / (W * H)
    C = K_I * K_ox * EPS0 * L * (W / xo + H / Ls)
    return R, C


def current_density_MA_cm2(current_mA, W_um, H_um):
    """배선 단면의 전류 밀도 [MA/cm^2]."""
    return current_mA * 1e-3 / (W_um * H_um * 1e-8) / 1e6
