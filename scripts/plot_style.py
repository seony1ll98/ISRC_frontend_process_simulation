"""모든 모듈 그래프가 공유하는 색상과 스타일.

범주형 색은 고정 순서로 쓰고, 온도나 시간처럼 순서가 있는 값은 한 색상의 명도 단계로 표현한다.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # blue, orange, aqua, yellow
BLUE_RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
BLUE_RAMP_6 = ["#86b6ef", "#5598e7", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": TEXT_2, "axes.labelcolor": TEXT, "xtick.color": TEXT_2, "ytick.color": TEXT_2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 2.0, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "legend.frameon": False, "figure.dpi": 150,
})


def save(fig, fig_dir, name):
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / name, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {fig_dir.parent.name}/{fig_dir.name}/{name}")
