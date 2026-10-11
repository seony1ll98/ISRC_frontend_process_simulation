"""포토리소그래피 모델.

- Rayleigh 해상도와 초점 심도(DOF)
- 근접(proximity) 노광: Fresnel 회절 (해석해와 각스펙트럼 수치 전파)
- 축소 투영(projection) 노광: Abbe 부분 결맞음 결상, 초점 이탈(defocus)
- 레지스트: contrast(gamma), CMTF, 문턱값(threshold) 모델로 CD 계산

공통 단위: 길이 [um], 파장 [um], 초점 이탈 [um]
"""

import numpy as np
from scipy.special import fresnel

WAVELENGTH_UM = {"g-line": 0.436, "i-line": 0.365, "KrF": 0.248, "ArF": 0.193}


# ---------------------------------------------------------------------------
# Rayleigh 해상도와 DOF (Plummer Ch. 5, Eq. 4, 5)
# ---------------------------------------------------------------------------
def rayleigh_resolution(wavelength_um, NA, k1=0.61):
    """해상도 R = k1 * lambda / NA [um]."""
    return k1 * wavelength_um / NA


def depth_of_focus(wavelength_um, NA, k2=0.5):
    """초점 심도 DOF = +/- k2 * lambda / NA^2 [um]. 반환값은 한쪽 폭이다."""
    return k2 * wavelength_um / NA**2


def k1_factor(half_pitch_um, wavelength_um, NA):
    """패턴의 공정 난이도 지표 k1 = half pitch * NA / lambda."""
    return half_pitch_um * NA / wavelength_um


# ---------------------------------------------------------------------------
# 레지스트 (Plummer Ch. 5, Eq. 11, 12)
# ---------------------------------------------------------------------------
def resist_gamma(D0, Df):
    """레지스트 contrast gamma = 1 / log10(Df / D0)."""
    return 1.0 / np.log10(Df / D0)


def resist_cmtf(gamma):
    """레지스트의 임계 MTF. 공간상의 contrast가 이보다 커야 패턴이 형성된다."""
    r = 10.0 ** (1.0 / gamma)
    return (r - 1.0) / (r + 1.0)


def image_contrast(I):
    """MTF 정의 (Imax - Imin) / (Imax + Imin)."""
    I = np.asarray(I, dtype=float)
    return float((I.max() - I.min()) / (I.max() + I.min()))


# ---------------------------------------------------------------------------
# 근접 노광 (Fresnel 회절)
# ---------------------------------------------------------------------------
def proximity_min_feature(wavelength_um, gap_um):
    """근접 노광의 최소 선폭 근사 W_min = sqrt(lambda * g) (Plummer Ch. 5, Eq. 10)."""
    return np.sqrt(wavelength_um * gap_um)


def fresnel_number(width_um, wavelength_um, gap_um):
    """N_F = W^2 / (lambda g). 1보다 충분히 크면 개구 모양이 그대로 전사된다."""
    return width_um**2 / (wavelength_um * gap_um)


def slit_fresnel_intensity(x_um, width_um, gap_um, wavelength_um):
    """폭 W 슬릿을 평면파로 비췄을 때 간격 g에서의 세기 (개구 없는 경우 = 1).

    1D Fresnel 회절 적분을 Fresnel 적분 C, S로 계산한다.
    """
    x = np.asarray(x_um, dtype=float)
    if gap_um == 0:
        return (np.abs(x) <= width_um / 2).astype(float)
    a = np.sqrt(2.0 / (wavelength_um * gap_um))
    S1, C1 = fresnel(a * (-width_um / 2 - x))
    S2, C2 = fresnel(a * (width_um / 2 - x))
    return 0.5 * ((C2 - C1) ** 2 + (S2 - S1) ** 2)


def angular_spectrum_propagate(x_um, field, gap_um, wavelength_um, paraxial=False):
    """1D 스칼라 장을 각스펙트럼 방법(FFT)으로 g만큼 전파한다.

    paraxial=False: 정확한 전달함수 exp(i 2 pi g sqrt(1/lambda^2 - f^2)), 소멸파는 감쇠
    paraxial=True : Fresnel(근축) 근사 exp(-i pi lambda g f^2). Fresnel 적분 해와 같은 근사다.
    """
    x = np.asarray(x_um, dtype=float)
    dx = x[1] - x[0]
    f = np.fft.fftfreq(x.size, d=dx)
    if paraxial:
        H = np.exp(-1j * np.pi * wavelength_um * gap_um * f**2)
    else:
        kz = np.sqrt((1.0 / wavelength_um) ** 2 - f.astype(complex) ** 2)
        H = np.exp(2j * np.pi * gap_um * (kz - 1.0 / wavelength_um))   # 공통 위상은 제거
    return np.fft.ifft(np.fft.fft(field) * H)


def slit_field(x_um, width_um):
    """폭 W 슬릿의 투과 진폭. 격자에 걸친 경계 칸은 걸친 비율만큼 투과시켜 이산화 오차를 줄인다."""
    x = np.asarray(x_um, dtype=float)
    dx = x[1] - x[0]
    return np.clip((width_um / 2 - np.abs(x)) / dx + 0.5, 0.0, 1.0).astype(complex)


def threshold_width(x_um, I, threshold, inside=True):
    """세기가 threshold 이상(inside=True)인 중앙 영역의 폭 [um]. 경계는 선형 보간한다."""
    x = np.asarray(x_um, dtype=float)
    I = np.asarray(I, dtype=float)
    above = I >= threshold if inside else I < threshold
    c = int(np.argmin(np.abs(x)))
    if not above[c]:
        return 0.0
    i = c
    while i > 0 and above[i - 1]:
        i -= 1
    j = c
    while j < x.size - 1 and above[j + 1]:
        j += 1
    if i == 0 or j == x.size - 1:
        return float("nan")

    def cross(k0, k1):
        return x[k0] + (threshold - I[k0]) * (x[k1] - x[k0]) / (I[k1] - I[k0])

    return float(cross(j, j + 1) - cross(i - 1, i))


# ---------------------------------------------------------------------------
# 축소 투영 노광 (Abbe 부분 결맞음 결상)
# ---------------------------------------------------------------------------
def line_space_orders(pitch_um, line_um, n_max):
    """선(불투명)이 x = 0 중앙에 있는 binary 마스크의 회절 차수 진폭 a_n.

    t(x) = 1 - rect(x / L) (주기 p)
    """
    n = np.arange(-n_max, n_max + 1)
    duty = line_um / pitch_um
    a = np.where(n == 0, 1.0 - duty, -np.sin(np.pi * n * duty) / (np.pi * np.where(n == 0, 1, n)))
    return n, a


def source_points(sigma, n_across=41):
    """원형 조명(반지름 sigma, NA 단위)을 균일 격자로 샘플링한다. sigma = 0이면 축상 한 점."""
    if sigma == 0:
        return np.zeros((1, 2))
    g = np.linspace(-sigma, sigma, n_across)
    sx, sy = np.meshgrid(g, g)
    keep = sx**2 + sy**2 <= sigma**2 + 1e-12
    return np.column_stack([sx[keep], sy[keep]])


def aerial_image(x_um, pitch_um, line_um, wavelength_um, NA, sigma=0.0, defocus_um=0.0, n_across=41):
    """주기적 line/space 패턴의 공간상(aerial image). 개구 없는 마스크의 세기 = 1.

    각 조명점에서 동공(|rho| <= 1)을 통과한 회절 차수만 간섭시키고(Abbe 방법),
    조명점마다 얻은 세기를 평균한다. 초점 이탈은 동공 위상
    phi = 2 pi z / lambda * (sqrt(1 - NA^2 rho^2) - 1)로 넣는다.
    """
    x = np.asarray(x_um, dtype=float)
    f0 = wavelength_um / (pitch_um * NA)            # 1차 회절의 동공 좌표 간격
    n_max = int(np.floor((1.0 + sigma) / f0)) + 1
    n, a = line_space_orders(pitch_um, line_um, n_max)
    src = source_points(sigma, n_across)
    phase_x = np.exp(2j * np.pi * np.outer(n, x) / pitch_um)          # (orders, x)
    I = np.zeros_like(x)
    for sx, sy in src:
        rho2 = (sx + n * f0) ** 2 + sy**2
        passed = rho2 <= 1.0 + 1e-12
        phi = 2 * np.pi * defocus_um / wavelength_um * (np.sqrt(np.clip(1 - NA**2 * rho2, 0, None)) - 1.0)
        coeff = np.where(passed, a * np.exp(1j * phi), 0.0)
        I += np.abs(coeff @ phase_x) ** 2
    return I / len(src)


def cutoff_pitch(wavelength_um, NA, sigma):
    """부분 결맞음 조명에서 1차 회절광을 받을 수 있는 최소 피치 lambda / (NA (1 + sigma))."""
    return wavelength_um / (NA * (1.0 + sigma))


def _circle_overlap_area(r1, r2, d):
    """반지름 r1, r2인 두 원의 중심 거리가 d일 때 겹치는 면적."""
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return np.pi * min(r1, r2) ** 2
    a1 = r1**2 * np.arccos((d**2 + r1**2 - r2**2) / (2 * d * r1))
    a2 = r2**2 * np.arccos((d**2 + r2**2 - r1**2) / (2 * d * r2))
    a3 = 0.5 * np.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    return a1 + a2 - a3


def dense_ls_contrast_analytic(pitch_um, line_um, wavelength_um, NA, sigma):
    """0차와 1차만 결상에 참여하는 영역의 contrast 해석해 (검증용).

    조건: sigma <= 1, 그리고 1 < f0 (= lambda / (p NA))  -> +1차와 -1차가 동시에 통과하지 못한다.
    조명 원(반지름 sigma) 중 +1차가 통과하는 비율 F를 원의 겹침 면적으로 계산하면
    I(x) = a0^2 + 2 a1^2 F + 4 a0 a1 F cos(2 pi x / p) 이다.
    """
    f0 = wavelength_um / (pitch_um * NA)
    if not (sigma <= 1.0 and f0 > 1.0):
        raise ValueError("해석해는 sigma <= 1, lambda/(p NA) > 1 영역에서만 성립합니다.")
    _, a = line_space_orders(pitch_um, line_um, 1)
    a_m1, a0, a1 = a
    F = _circle_overlap_area(sigma, 1.0, f0) / (np.pi * sigma**2)
    dc = a0**2 + 2 * a1**2 * F
    amp = abs(4 * a0 * a1 * F)
    return amp / dc


def coherent_defocus_contrast_analytic(pitch_um, line_um, wavelength_um, NA, defocus_um):
    """결맞음 조명(sigma = 0), 0차와 +/-1차만 통과할 때 초점 이탈에 따른 contrast 해석해 (검증용).

    I(x) = a0^2 + 2 a1^2 + 4 a0 a1 cos(dphi) cos(kx) + 2 a1^2 cos(2kx)
    """
    f0 = wavelength_um / (pitch_um * NA)
    if not (f0 <= 1.0 < 2 * f0):
        raise ValueError("0차와 +/-1차만 통과하는 피치에서만 성립합니다.")
    _, a = line_space_orders(pitch_um, line_um, 1)
    a0, a1 = a[1], a[2]
    dphi = 2 * np.pi * defocus_um / wavelength_um * (np.sqrt(1 - (wavelength_um / pitch_um) ** 2) - 1.0)
    x = np.linspace(-pitch_um / 2, pitch_um / 2, 4001)
    k = 2 * np.pi / pitch_um
    I = a0**2 + 2 * a1**2 + 4 * a0 * a1 * np.cos(dphi) * np.cos(k * x) + 2 * a1**2 * np.cos(2 * k * x)
    return image_contrast(I)


# ---------------------------------------------------------------------------
# 레지스트 문턱값 모델과 공정 창
# ---------------------------------------------------------------------------
def resist_line_cd(x_um, I, dose, threshold=1.0):
    """양성 레지스트에서 남는 선(x = 0 중앙)의 폭 [um].

    dose * I(x) < threshold 인 곳은 노광이 부족해 레지스트가 남는다.
    x는 한 주기 [-p/2, p/2]를 덮어야 한다.
    """
    x = np.asarray(x_um, dtype=float)
    E = dose * np.asarray(I, dtype=float)
    remain = E < threshold
    c = int(np.argmin(np.abs(x)))
    if not remain[c]:
        return 0.0
    if remain.all():
        return float(x[-1] - x[0])
    j = c
    while j < x.size - 1 and remain[j + 1]:
        j += 1
    if j == x.size - 1:
        return float(x[-1] - x[0])
    xr = x[j] + (threshold - E[j]) * (x[j + 1] - x[j]) / (E[j + 1] - E[j])
    return float(2 * xr)                             # 대칭 패턴


def dose_to_size(x_um, I, target_cd_um, threshold=1.0):
    """목표 CD가 나오는 노광량(문턱값 대비 배수)을 구한다."""
    from scipy.optimize import brentq
    f = lambda d: resist_line_cd(x_um, I, d, threshold) - target_cd_um
    lo, hi = threshold / I.max() * 1.0001, threshold / I.min() * 0.9999
    return brentq(f, lo, hi)
