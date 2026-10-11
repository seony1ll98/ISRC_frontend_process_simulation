"""박막 증착 모델.

- CVD 성장 속도: Grove 모델(표면 반응과 물질 전달의 직렬 결합)
- 트렌치 안의 단차 피복(step coverage): 입사각 분포 cos^n과 부착 확률 S_c, 재방출(re-emission)

공통 단위: 길이 [um], 온도 [C]
"""

import numpy as np

from .constants import K_B_EV, celsius_to_kelvin


# ---------------------------------------------------------------------------
# CVD: Grove 모델 (Plummer Ch. 9, Eq. 8~12)
# ---------------------------------------------------------------------------
def grove_growth_rate(ks, hG, CT_over_N=1.0, Y=1.0):
    """v = (ks hG / (ks + hG)) * (C_T / N) * Y. ks와 hG는 같은 단위(cm/s)."""
    ks = np.asarray(ks, dtype=float)
    return ks * hG / (ks + hG) * CT_over_N * Y


def surface_rate_constant(T_C, Ea_eV, T_transition_C, hG=1.0):
    """ks = k0 exp(-Ea/kT). k0는 T_transition에서 ks = hG가 되도록 정한다(물질 전달과 반응 속도가 같은 온도)."""
    kT = K_B_EV * celsius_to_kelvin(np.asarray(T_C, dtype=float))
    kT0 = K_B_EV * celsius_to_kelvin(T_transition_C)
    return hG * np.exp(-Ea_eV / kT + Ea_eV / kT0)


def transition_temperature_shift(T_transition_C, Ea_eV, hG_ratio):
    """hG가 hG_ratio배가 될 때 전이 온도(ks = hG)가 어디로 이동하는지 [C]."""
    inv_T = 1.0 / celsius_to_kelvin(T_transition_C) - K_B_EV * np.log(hG_ratio) / Ea_eV
    return 1.0 / inv_T - 273.15


def temperature_sensitivity(T_C, Ea_eV):
    """반응 율속 영역에서 성장 속도의 온도 민감도 d(ln v)/dT = Ea / kT^2 [1/K]."""
    T = celsius_to_kelvin(T_C)
    return Ea_eV / (K_B_EV * T**2)


def lpcvd_hG_gain(pressure_torr, boundary_layer_growth=7.0):
    """압력을 낮출 때 물질 전달 계수의 증가 배수. D_G ∝ 1/P, 경계층 두께는 boundary_layer_growth배."""
    return (760.0 / pressure_torr) / boundary_layer_growth


# ---------------------------------------------------------------------------
# 트렌치 단차 피복
# ---------------------------------------------------------------------------
def trench_elements(width, depth, n_per_side=40, n_bottom=None):
    """2D 트렌치(깊이 방향으로 무한히 긴 홈) 내부 표면을 선분 요소로 나눈다.

    좌표: 트렌치 입구가 y = 0, 바닥이 y = -depth, x는 -w/2 ~ w/2.
    반환: 요소 양 끝점 (p0, p1), 중점, 길이, 바깥(빈 공간 쪽) 법선, 요소 종류
    """
    if n_bottom is None:
        n_bottom = max(int(round(n_per_side * width / depth)), 4)
    w2 = width / 2
    segs, kinds = [], []
    ys = np.linspace(0, -depth, n_per_side + 1)
    for a, b in zip(ys[:-1], ys[1:]):                       # 왼쪽 벽, 위에서 아래로
        segs.append(((-w2, a), (-w2, b)))
        kinds.append("left")
    xs = np.linspace(-w2, w2, n_bottom + 1)
    for a, b in zip(xs[:-1], xs[1:]):                       # 바닥, 왼쪽에서 오른쪽으로
        segs.append(((a, -depth), (b, -depth)))
        kinds.append("bottom")
    for a, b in zip(ys[::-1][:-1], ys[::-1][1:]):           # 오른쪽 벽, 아래에서 위로
        segs.append(((w2, a), (w2, b)))
        kinds.append("right")
    segs = np.array(segs, dtype=float)
    p0, p1 = segs[:, 0], segs[:, 1]
    mid = 0.5 * (p0 + p1)
    length = np.linalg.norm(p1 - p0, axis=1)
    normal = np.zeros_like(mid)
    kinds = np.array(kinds)
    normal[kinds == "left"] = (1, 0)
    normal[kinds == "right"] = (-1, 0)
    normal[kinds == "bottom"] = (0, 1)
    return {"p0": p0, "p1": p1, "mid": mid, "length": length, "normal": normal, "kind": kinds}


def view_factor_segments(a0, a1, b0, b1):
    """2D 선분 A에서 선분 B로의 view factor (crossed-string 방법, 가림이 없을 때 정확)."""
    d = lambda p, q: np.linalg.norm(np.asarray(p) - np.asarray(q))
    # 끝점 순서에 따라 '교차'와 '비교차' 연결이 바뀌므로 절댓값을 쓴다.
    return abs((d(a0, b1) + d(a1, b0)) - (d(a0, b0) + d(a1, b1))) / (2 * d(a0, a1))


def view_factor_matrix(el):
    """요소 i에서 요소 k로의 view factor 행렬. 직사각형 트렌치 내부는 서로 가리지 않는다."""
    n = len(el["length"])
    V = np.zeros((n, n))
    for i in range(n):
        for k in range(n):
            if i == k or el["kind"][i] == el["kind"][k]:
                continue                                    # 같은 평면 위 요소끼리는 서로 보지 않는다
            V[i, k] = view_factor_segments(el["p0"][i], el["p1"][i], el["p0"][k], el["p1"][k])
    return V


def view_factor_to_opening(el, width):
    """각 요소에서 트렌치 입구(y = 0)로의 view factor."""
    w2 = width / 2
    return np.array([view_factor_segments(a, b, (w2, 0.0), (-w2, 0.0)) for a, b in zip(el["p0"], el["p1"])])


def direct_flux(el, width, n, n_theta=300, n_phi=240):
    """입구 위에서 cos^n 분포로 들어오는 입자가 각 요소(중점)에 직접 도달하는 플럭스 (평평한 면 = 1).

    트렌치가 깊이 방향으로 무한히 길어, 입자 방향을 단면에 투영한 각도만으로 가림 여부가 정해진다.
    """
    th = (np.arange(n_theta) + 0.5) * (np.pi / 2) / n_theta
    ph = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    dOm = np.sin(TH) * (np.pi / 2 / n_theta) * (2 * np.pi / n_phi)
    # u: 입자가 온 방향(위쪽). 단면 성분 (ux, uz), 깊이 방향 uy
    ux, uz = np.sin(TH) * np.cos(PH), np.cos(TH)
    radiance = (n + 1) / (2 * np.pi) * np.cos(TH) ** (n - 1)  # 수평면 플럭스 = 1이 되도록 정규화
    w2 = width / 2
    out = np.zeros(len(el["length"]))
    for i, (p, m) in enumerate(zip(el["mid"], el["normal"])):
        cos_in = ux * m[0] + uz * m[1]
        # 투영 방향의 직선이 y = 0에 닿는 x 위치 (uz > 0이므로 항상 위로 올라간다)
        x_hit = p[0] + ux / uz * (0.0 - p[1])
        visible = (x_hit > -w2) & (x_hit < w2)
        out[i] = np.sum(radiance * np.clip(cos_in, 0, None) * visible * dOm)
    return out


def solve_trench_flux(el, width, n=1, Sc=1.0, V=None):
    """재방출을 포함한 각 요소의 전체 도달 플럭스 F (평평한 면의 입사 플럭스 = 1).

    F_i = F_direct,i + sum_k V_ik (1 - S_c) F_k   (Plummer Ch. 9, Eq. 19, 23)
    재방출은 cos 분포(입사 방향을 기억하지 않음)로 둔다. 증착 속도는 S_c F에 비례한다.
    """
    if V is None:
        V = view_factor_matrix(el)
    Fd = direct_flux(el, width, n)
    A = np.eye(len(Fd)) - (1.0 - Sc) * V
    F = np.linalg.solve(A, Fd)
    return F, Fd, V


def step_coverage_metrics(el, F, Sc=None):
    """증착 속도(S_c F)를 평평한 윗면(S_c x 1)으로 나눈 상대 두께 지표. S_c는 약분되어 F만 쓴다."""
    rate = F                                                # S_c는 분자와 분모에 같이 들어가 약분된다
    kind, mid = el["kind"], el["mid"]
    bottom = rate[kind == "bottom"]
    left = rate[kind == "left"]
    ly = mid[kind == "left", 1]
    order = np.argsort(-ly)                                 # 위에서 아래로
    left = left[order]
    return {
        "bottom_center": float(bottom[len(bottom) // 2]),
        "bottom_corner": float(bottom[0]),
        "sidewall_top": float(left[0]),
        "sidewall_mid": float(left[len(left) // 2]),
        "sidewall_bottom": float(left[-1]),
        "min_over_trench": float(rate.min()),
    }
