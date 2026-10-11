"""건식 식각 모델.

- 식각 지표: 이방성(anisotropy), 선택비(selectivity), 과식각(over-etch)과 하부막 손실
- 식각 속도 모델: 선형 모델, 포화 흡착(이온 보조) 모델 (Plummer Ch. 10, Eq. 7, 8)
- 이온의 입사각 분포 cos^n(theta)에 따른 기울어진 면의 이온 플럭스
- 2D 단면 형상 진화: level set 방법

공통 단위: 길이 [um], 시간 [min], 식각 속도 [um/min]
"""

import numpy as np
import skfmm
from scipy.ndimage import distance_transform_edt


# ---------------------------------------------------------------------------
# 식각 지표 (Plummer Ch. 10, Eq. 3)
# ---------------------------------------------------------------------------
def anisotropy(r_lateral, r_vertical):
    """A_f = 1 - r_lat / r_ver. 등방성이면 0, 완전 이방성이면 1."""
    return 1.0 - r_lateral / r_vertical


def anisotropy_from_bias(bias_um, depth_um):
    """마스크 아래로 파고든 폭(bias, 한쪽) b와 식각 깊이 d로 A_f = 1 - b/d."""
    return 1.0 - bias_um / depth_um


def selectivity(r_film, r_other):
    """선택비 S = r_film / r_other."""
    return r_film / r_other


# ---------------------------------------------------------------------------
# 과식각과 선택비 설계
# ---------------------------------------------------------------------------
def clear_time_window(thickness_um, rate_um_min, thickness_nonuniformity, rate_nonuniformity):
    """웨이퍼 안에서 막이 가장 먼저 뚫리는 시간과 가장 늦게 뚫리는 시간 [min].

    두께와 식각 속도가 각각 +/- 비율만큼 퍼져 있다고 본다(최악 조합).
    """
    t_first = thickness_um * (1 - thickness_nonuniformity) / (rate_um_min * (1 + rate_nonuniformity))
    t_last = thickness_um * (1 + thickness_nonuniformity) / (rate_um_min * (1 - rate_nonuniformity))
    return t_first, t_last


def underlayer_loss(thickness_um, rate_um_min, selectivity_to_under, thickness_nonuniformity,
                    rate_nonuniformity, overetch_fraction):
    """최악 위치(가장 먼저 뚫린 곳)에서 하부막이 깎이는 두께 [um].

    공정 시간은 '가장 늦게 뚫리는 시간 x (1 + 과식각 비율)'로 둔다.
    가장 먼저 뚫린 곳은 그 사이 하부막이 노출되어 (식각 속도 / 선택비)로 깎인다.
    """
    t_first, t_last = clear_time_window(thickness_um, rate_um_min, thickness_nonuniformity, rate_nonuniformity)
    t_total = t_last * (1 + overetch_fraction)
    r_fast = rate_um_min * (1 + rate_nonuniformity)
    return (t_total - t_first) * r_fast / selectivity_to_under


def required_selectivity(thickness_um, max_loss_um, thickness_nonuniformity, rate_nonuniformity, overetch_fraction):
    """하부막 손실을 max_loss 이하로 유지하는 최소 선택비. 식각 속도에 무관하다."""
    loss_at_S1 = underlayer_loss(thickness_um, 1.0, 1.0, thickness_nonuniformity, rate_nonuniformity,
                                 overetch_fraction)
    return loss_at_S1 / max_loss_um


def required_mask_thickness(thickness_um, selectivity_to_mask, thickness_nonuniformity, rate_nonuniformity,
                            overetch_fraction, remaining_um=0.0):
    """공정이 끝날 때까지 마스크(PR)가 remaining 이상 남는 데 필요한 초기 마스크 두께 [um]."""
    _, t_last = clear_time_window(thickness_um, 1.0, thickness_nonuniformity, rate_nonuniformity)
    t_total = t_last * (1 + overetch_fraction)
    return t_total * (1 + rate_nonuniformity) / selectivity_to_mask + remaining_um


# ---------------------------------------------------------------------------
# 식각 속도 모델 (Plummer Ch. 10, Eq. 7, 8)
# ---------------------------------------------------------------------------
def rate_linear(chem, ion):
    """선형 모델: 화학 성분(S_c K_f F_c / N)과 이온 성분(K_i F_i / N)이 독립적으로 더해진다."""
    return np.asarray(chem, dtype=float) + np.asarray(ion, dtype=float)


def rate_saturation(chem, ion):
    """포화 흡착 모델: 두 과정이 직렬로 일어나 느린 쪽이 속도를 제한한다.

    ER = 1 / (1/ion + 1/chem). 어느 한쪽이 0이면 식각되지 않는다.
    """
    chem = np.asarray(chem, dtype=float)
    ion = np.asarray(ion, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = chem * ion / (chem + ion)
    return np.where((chem > 0) & (ion > 0), out, 0.0)


# ---------------------------------------------------------------------------
# 이온 입사각 분포와 기울어진 면의 플럭스
# ---------------------------------------------------------------------------
def directional_flux_factor(tilt_rad, n, n_theta=400, n_phi=200):
    """입사각 분포가 cos^n(theta)인 입자가 기울어진 면에 주는 플럭스 (수평면 = 1).

    cos^n(theta)는 수평면에 도달하는 단위 입체각당 플럭스의 분포다(Plummer Ch. 10).
    면의 법선이 연직(아래 방향)에서 tilt만큼 기울어져 있을 때, 그 면이 받는 플럭스를
    입체각 적분으로 계산한다. n = 1이면 해석해 (1 + cos tilt) / 2, n이 크면 cos tilt에 접근한다.
    """
    tilt = np.atleast_1d(np.asarray(tilt_rad, dtype=float))
    th = (np.arange(n_theta) + 0.5) * (np.pi / 2) / n_theta
    ph = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    dOmega = np.sin(TH) * (np.pi / 2 / n_theta) * (2 * np.pi / n_phi)
    radiance = np.cos(TH) ** (n - 1)                     # 수평면 플럭스 cos^n = 복사휘도 x cos
    # 입자 진행 방향 w (z축 = 아래 방향 기준), 면 법선 m = (sin tilt, 0, cos tilt)
    wx, wz = np.sin(TH) * np.cos(PH), np.cos(TH)
    flat = np.sum(radiance * wz * dOmega)
    out = np.array([np.sum(radiance * np.clip(wx * np.sin(a) + wz * np.cos(a), 0, None) * dOmega) for a in tilt])
    out = out / flat
    return out if np.ndim(tilt_rad) else float(out[0])


# ---------------------------------------------------------------------------
# 2D 단면 형상 진화 (level set)
# ---------------------------------------------------------------------------
VOID, MASK, FILM, UNDER = 0, 1, 2, 3


def make_structure(opening_um, film_um, mask_um=0.3, under_um=0.3, width_um=2.0, top_um=0.2, dx_um=0.005):
    """마스크(PR) / 식각할 막 / 하부막이 쌓인 구조를 격자로 만든다. 막 윗면이 y = 0이다."""
    x = np.arange(-width_um / 2, width_um / 2 + dx_um / 2, dx_um)
    y = np.arange(-(film_um + under_um), mask_um + top_um + dx_um / 2, dx_um)
    X, Y = np.meshgrid(x, y, indexing="xy")
    mat = np.full(X.shape, VOID, dtype=np.int8)
    mat[(Y >= 0) & (Y < mask_um) & (np.abs(X) >= opening_um / 2)] = MASK
    mat[(Y < 0) & (Y >= -film_um)] = FILM
    mat[Y < -film_um] = UNDER
    return x, y, mat


def _signed_distance(void, dx):
    """빈 공간은 음수, 고체는 양수인 부호 거리."""
    return (distance_transform_edt(~void) - distance_transform_edt(void)) * dx + np.where(void, 0.5, -0.5) * dx


def _reinitialize(phi, dx):
    """경계(phi = 0)의 위치를 유지한 채 phi를 부호 거리로 되돌린다(fast marching, 2차 정확도)."""
    return skfmm.distance(phi, dx=float(dx), order=2)


def simulate_profile(x, y, mat, R_chem, R_ion, t_end, model="linear", ion_n=None,
                     sel_under=np.inf, sel_mask=np.inf, cfl=0.45, n_snapshots=5):
    """식각 단면 형상을 시간에 따라 진화시킨다.

    각 표면 지점의 법선 속도 V는 법선 방향에 따라 정한다.
      - 화학 성분 R_chem: 반응성 중성종은 부착 확률이 낮아(S_c ~ 0.01) 여러 번 튕기며
        모든 방향의 면에 고르게 도달한다고 보고 방향과 무관하게 둔다.
      - 이온 성분 R_ion: 위에서 아래로 입사한다. ion_n이 None이면 완전히 연직,
        숫자면 cos^n 분포로 directional_flux_factor를 곱한다.
      - model: "linear"는 두 성분의 합, "saturation"은 직렬 결합.
    하부막은 속도를 sel_under로, 마스크는 sel_mask로 나눈다(inf면 식각되지 않음).
    가림(shadowing), 재증착, 대전 효과는 다루지 않는다.

    Returns: 스냅숏 시간 목록과 각 시점의 빈 공간(bool 배열) 목록
    """
    dx = x[1] - x[0]
    void = mat == VOID
    phi = _signed_distance(void, dx)
    if ion_n is not None:
        a_tab = np.linspace(0, np.pi, 721)
        g_tab = directional_flux_factor(a_tab, ion_n)

    def speed(cos_tilt):
        """법선과 연직 아래 방향 사이각의 cos에 따른 법선 속도."""
        if ion_n is None:
            g = np.clip(cos_tilt, 0, None)
        else:
            g = np.interp(np.arccos(np.clip(cos_tilt, -1, 1)), a_tab, g_tab)
        ion = R_ion * g
        return rate_saturation(R_chem, ion) if model == "saturation" else rate_linear(R_chem, ion)

    def hamiltonian(px, py):
        """H(p) = |p| V(p/|p|). 법선 p/|p|는 고체 쪽을 향한다."""
        norm = np.hypot(px, py)
        return norm * speed(-py / np.maximum(norm, 1e-12))

    # Lax-Friedrichs 소산 계수: H는 1차 동차이므로 단위원 위에서 |dH/dp|의 최댓값을 구하면 된다.
    th = np.linspace(0, 2 * np.pi, 2881)
    ux, uy, h = np.cos(th), np.sin(th), 1e-4
    a_x = np.max(np.abs(hamiltonian(ux + h, uy) - hamiltonian(ux - h, uy)) / (2 * h))
    a_y = np.max(np.abs(hamiltonian(ux, uy + h) - hamiltonian(ux, uy - h)) / (2 * h))
    a_x, a_y = 1.05 * a_x + 1e-12, 1.05 * a_y + 1e-12

    scale = np.ones_like(phi)
    scale[mat == UNDER] = 1.0 / sel_under
    scale[mat == MASK] = 1.0 / sel_mask

    def extend_scale(phi):
        """빈 공간 격자점에는 가장 가까운 고체 격자점의 재료 속도 배율을 준다.

        경계 위치는 양쪽 격자점의 phi로 정해지므로, 빈 공간 쪽 격자점도 맞닿은 재료의 속도로
        움직여야 마스크처럼 식각되지 않는 재료가 경계에서 깎이지 않는다.
        """
        solid = phi >= 0
        _, (iy, ix) = distance_transform_edt(~solid, return_indices=True)
        return scale[iy, ix]

    ext = extend_scale(phi)
    dt = cfl * dx / (a_x + a_y)
    n_steps = int(np.ceil(t_end / dt))
    dt = t_end / n_steps
    snap_steps = set(np.linspace(0, n_steps, n_snapshots).round().astype(int))
    times, snaps = [], []
    for k in range(n_steps + 1):
        if k in snap_steps:
            times.append(k * dt)
            snaps.append(phi < 0)
        if k == n_steps:
            break
        p = np.pad(phi, 1, mode="edge")
        dxm = (p[1:-1, 1:-1] - p[1:-1, :-2]) / dx
        dxp = (p[1:-1, 2:] - p[1:-1, 1:-1]) / dx
        dym = (p[1:-1, 1:-1] - p[:-2, 1:-1]) / dx
        dyp = (p[2:, 1:-1] - p[1:-1, 1:-1]) / dx
        # phi_t + H(grad phi) = 0, 국소 Lax-Friedrichs(LLF)
        # 소산 계수를 격자점마다 한쪽 차분 조합에서의 |dH/dp| 최댓값으로 정해, 곧은 측벽처럼
        # 실제 속도가 작은 곳에 불필요한 수치 확산이 생기지 않게 한다.
        H = hamiltonian(0.5 * (dxm + dxp), 0.5 * (dym + dyp))
        lx = np.zeros_like(phi)
        ly = np.zeros_like(phi)
        for qx in (dxm, dxp):
            for qy in (dym, dyp):
                nrm = np.hypot(qx, qy) + 1e-12
                ux, uy = qx / nrm, qy / nrm
                lx = np.maximum(lx, np.abs(hamiltonian(ux + 1e-4, uy) - hamiltonian(ux - 1e-4, uy)) / 2e-4)
                ly = np.maximum(ly, np.abs(hamiltonian(ux, uy + 1e-4) - hamiltonian(ux, uy - 1e-4)) / 2e-4)
        diss = 0.5 * lx * (dxp - dxm) + 0.5 * ly * (dyp - dym)
        phi = phi - dt * ext * (H - diss)
        if k % 5 == 4:
            phi = _reinitialize(phi, dx)
            ext = extend_scale(phi)
    return times, snaps


def profile_metrics(x, y, void, at_depth_um=0.1):
    """식각 깊이(개구 중앙)와, 막 윗면에서 at_depth만큼 아래 위치의 식각 폭 절반을 계산한다.

    bias는 마스크 모서리 바로 아래의 특이점을 피해 막 윗면보다 조금 아래에서 잰다.
    """
    dx = x[1] - x[0]
    ix0 = int(np.argmin(np.abs(x)))
    col = void[:, ix0]
    below = np.nonzero(col & (y < 0))[0]
    depth = float(-(y[below.min()] - dx / 2)) if below.size else 0.0
    iy = int(np.argmin(np.abs(y + at_depth_um)))
    xs = x[void[iy] & (x >= 0)]
    half_width = float(xs.max() + dx / 2) if xs.size else 0.0
    return {"depth_um": depth, "half_width_um": half_width}
