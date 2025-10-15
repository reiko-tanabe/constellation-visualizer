# plot_constellation_plotly.py
# 使い方:
#   python plot_constellation_plotly.py stars.csv edges.json out.html [--no-labels] [--title=...]
#   python plot_constellation_plotly.py stars.csv out.html [--no-labels] [--title=...]

import sys, csv, json, math, os, unicodedata, re
from typing import List, Tuple, Dict, Any, Optional
import plotly.graph_objects as go
import numpy as np  # サイズ計算用

# ---------- 角度/座標ユーティリティ ----------
def ra_hms_to_deg(h: float, m: float, s: float) -> float:
    return (h + m/60.0 + s/3600.0) * 15.0

def dec_dms_to_deg(sign: int, d: float, m: float, s: float) -> float:
    val = abs(d) + m/60.0 + s/3600.0
    return val if sign >= 0 else -val

def sph_to_cart(ra_deg: float, dec_deg: float, dist: float):
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    x = dist * math.cos(dec) * math.cos(ra)
    y = dist * math.cos(dec) * math.sin(ra)
    z = dist * math.sin(dec)
    return x, y, z

# ---------- 文字列正規化＆ベース名 ----------
def norm(s: str) -> str:
    if s is None:
        return ""
    return unicodedata.normalize("NFKC", s).strip()

def base_key(s: str) -> str:
    """'Hamal (α Ari)' → 'hamal' のように括弧以降を落として照合キー化"""
    s2 = norm(s)
    s2 = re.sub(r"\s*\(.*?\)\s*$", "", s2)  # 末尾の括弧表記を削除
    s2 = re.sub(r"\s+", " ", s2)            # 連続空白を1つに
    return s2.lower()

# ---------- 質量ざっくり推定（視等級+距離） ----------
# 1 pc ≈ 3.26 ly, 太陽の可視絶対等級 ~4.83
def estimate_mass_from_mag_dist(mag: Optional[float], dist_ly: Optional[float]) -> Optional[float]:
    if mag is None or dist_ly is None or dist_ly <= 0:
        return None
    try:
        dist_pc = dist_ly / 3.26
        M = mag - 5 * (math.log10(dist_pc) - 1.0)
        L = 10 ** ((4.83 - M) / 2.5)
        mass = L ** (1.0 / 3.5)  # 主系列近似
        if not np.isfinite(mass) or mass <= 0:
            return None
        return float(mass)
    except Exception:
        return None

# ---------- CSV/JSON ロード ----------
def _float_or_none(v: Any) -> Optional[float]:
    try:
        if v is None:
            return None
        s = str(v).strip()
        if s == "" or s.lower() in ("nan", "none"):
            return None
        return float(s)
    except Exception:
        return None

def load_stars_csv(csv_path: str) -> List[Dict[str, Any]]:
    """
    2系統のCSVを自動判別:
      A) 旧: ra_h, ra_m, ra_s, dec_sign, dec_d, dec_m, dec_s, dist_ly, mag ...
      B) 新: ra_deg, dec_deg, dist_ly, mag, (spec_type など)
    """
    stars = []
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = [c.strip() for c in (reader.fieldnames or [])]

        use_gaia_schema = ("ra_deg" in fieldnames and "dec_deg" in fieldnames and "dist_ly" in fieldnames)
        use_hms_schema  = ("ra_h"  in fieldnames and "ra_m"   in fieldnames and "ra_s"     in fieldnames and
                           "dec_sign" in fieldnames and "dec_d" in fieldnames and "dec_m" in fieldnames and "dec_s" in fieldnames)

        if not (use_gaia_schema or use_hms_schema):
            raise ValueError("CSVの列構成を判別できませんでした。ra_deg/dec_deg/dist_ly か、ra_h/.../dec_s を含めてください。")

        for row in reader:
            name = row.get("name", "").strip()
            mag  = _float_or_none(row.get("mag"))

            # ---- 角度と距離を取り出し ----
            if use_gaia_schema:
                ra_deg  = _float_or_none(row.get("ra_deg"))
                dec_deg = _float_or_none(row.get("dec_deg"))
                dist    = _float_or_none(row.get("dist_ly"))
            else:
                ra_deg  = ra_hms_to_deg(float(row["ra_h"]), float(row["ra_m"]), float(row["ra_s"]))
                dec_deg = dec_dms_to_deg(int(row["dec_sign"]), float(row["dec_d"]), float(row["dec_m"]), float(row["dec_s"]))
                dist    = _float_or_none(row.get("dist_ly"))

            # 距離が無いと3D座標が出せないのでスキップ（警告）
            if dist is None or dist <= 0:
                print(f"[WARN] skip (no dist_ly): {name}")
                continue

            # 任意列
            radius_sun = _float_or_none(row.get("radius_sun"))
            mass_sun   = _float_or_none(row.get("mass_sun"))
            spec_type  = row.get("spec_type", "").strip() if "spec_type" in row else ""

            # 位置計算
            x, y, z = sph_to_cart(ra_deg, dec_deg, dist)

            # 質量が未入力なら推定（ざっくり）
            if mass_sun is None:
                mass_sun = estimate_mass_from_mag_dist(mag, dist)

            stars.append({
                "name": name,
                "name_norm": norm(name),
                "ra_deg": ra_deg,
                "dec_deg": dec_deg,
                "dist_ly": dist,
                "mag": mag,
                "spec_type": spec_type or None,
                "radius_sun": radius_sun,
                "mass_sun": mass_sun,
                "x": x, "y": y, "z": z
            })
    return stars

def load_edges_json(edges_path: Optional[str]) -> List[Tuple[str, str]]:
    if not edges_path or not os.path.exists(edges_path):
        return []
    with open(edges_path, encoding="utf-8") as f:
        raw = json.load(f)
    edges = []
    for item in raw:
        if isinstance(item, (list, tuple)) and len(item) == 2:
            a, b = item
            edges.append((str(a), str(b)))
    return edges

# ---------- AU円（未使用なら呼ばれません） ----------
AU_LY = 1.58125074e-5  # 1 AU in light-years
def add_au_circle(fig, au_radius: float, exaggerate: float = 1.0, z_plane: float = 0.0, name: Optional[str] = None):
    r_ly = au_radius * AU_LY * exaggerate
    t = np.linspace(0, 2*np.pi, 240)
    x = r_ly * np.cos(t)
    y = r_ly * np.sin(t)
    z = np.full_like(x, z_plane)
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode="lines",
        line=dict(width=2, dash="dot"),
        name=name or (f"{au_radius} AU" + (f" ×{int(exaggerate)}" if exaggerate != 1.0 else "")),
        hoverinfo="skip"
    ))

# ---------- 3D図をHTMLへ ----------
def plot_constellation_html(csv_path: str, out_html: str, edges_path: Optional[str] = None,
                            title: Optional[str] = None, show_labels: bool = True):
    stars = load_stars_csv(csv_path)
    edges_raw = load_edges_json(edges_path)
    if not title:
        title = os.path.splitext(os.path.basename(csv_path))[0]

    # デバッグ出力
    print(f"[INFO] stars loaded: {len(stars)} from {csv_path}")
    if edges_path:
        print(f"[INFO] edges file: {edges_path}")
    print(f"[INFO] edges loaded: {len(edges_raw)}")

    # 名前→インデックス（正規化キー）
    name_to_index: Dict[str, int] = {s["name_norm"]: i for i, s in enumerate(stars)}

    # base_key（括弧落としキー）→ index（重複がなければ登録）
    base_to_index: Dict[str, int] = {}
    _dup: Dict[str, int] = {}
    for i, s in enumerate(stars):
        bk = base_key(s["name"])
        _dup[bk] = _dup.get(bk, 0) + 1
        if bk not in base_to_index:
            base_to_index[bk] = i
    # 重複するベース名は削除（曖昧回避）
    for bk, cnt in _dup.items():
        if cnt > 1 and bk in base_to_index:
            del base_to_index[bk]

    xs = [s["x"] for s in stars]
    ys = [s["y"] for s in stars]
    zs = [s["z"] for s in stars]
    texts = [s["name"] for s in stars]

    # ========== 太陽=1 の半径比でマーカーサイズ ==========
    def radius_ratio_vs_sun(star: Dict[str, Any]) -> Optional[float]:
        r = star.get("radius_sun")
        if r is not None and r > 0:
            return float(r)
        m = star.get("mass_sun")
        if m is not None and m > 0:
            return float(max(0.1, m ** 0.8))  # 最小0.1で暴れ抑制
        rr_mass = estimate_mass_from_mag_dist(star.get("mag"), star.get("dist_ly"))
        if rr_mass is not None and rr_mass > 0:
            return float(max(0.1, rr_mass ** 0.8))
        return None

    SUN_MARKER_PX = 5.0
    MAX_MARKER_PX = 48.0
    MIN_MARKER_PX = 3.0

    radius_ratios = [radius_ratio_vs_sun(s) for s in stars]

    def marker_px_from_ratio(ratio: Optional[float], fallback_mag: Optional[float]) -> float:
        if ratio is not None:
            return float(np.clip(SUN_MARKER_PX * ratio, MIN_MARKER_PX, MAX_MARKER_PX))
        if fallback_mag is not None:
            return float(np.clip(10.0 - (fallback_mag - 1.0), MIN_MARKER_PX, 18.0))
        return SUN_MARKER_PX

    sizes = [marker_px_from_ratio(rr, s.get("mag")) for rr, s in zip(radius_ratios, stars)]

    fig = go.Figure()

    # customdata: [mag, radius_sun, mass_sun, radius_ratio, spec_type]
    def _fmt(v):
        if v is None:
            return ""
        if isinstance(v, (int, float)):
            return f"{v:.3f}"
        return v
    customdata = [[_fmt(s.get("mag")), _fmt(s.get("radius_sun")), _fmt(s.get("mass_sun")),
                   _fmt(rr if rr is not None else None), _fmt(s.get("spec_type"))]
                  for s, rr in zip(stars, radius_ratios)]

    # 星（ラベルとホバー）
    mode = "markers+text" if show_labels else "markers"
    fig.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs,
        mode=mode,
        text=texts if show_labels else None,
        textposition="top center",
        marker=dict(size=sizes),
        hovertemplate=(
            "<b>%{text}</b>"
            "<br>x=%{x:.3f} ly, y=%{y:.3f} ly, z=%{z:.3f} ly"
            "<br>mag=%{customdata[0]}  radius☉=%{customdata[1]}  mass☉=%{customdata[2]}"
            "<br>radius ratio (vs Sun)=%{customdata[3]}  spec=%{customdata[4]}"
            "<extra></extra>"
        ),
        customdata=customdata,
        name="Stars"
    ))

    # 太陽（原点）
    fig.add_trace(go.Scatter3d(
        x=[0], y=[0], z=[0],
        mode="markers+text" if show_labels else "markers",
        text=["Sun (×1)"] if show_labels else None,
        textposition="bottom center",
        marker=dict(size=SUN_MARKER_PX, color="red"),
        name="Sun"
    ))

    # ------- 星座線（頑丈な名前解決） -------
    def resolve_index(name_raw: str) -> Optional[int]:
        n = norm(name_raw)
        if n in name_to_index:
            return name_to_index[n]
        bk = base_key(name_raw)
        if bk in base_to_index:
            return base_to_index[bk]
        return None

    missing_edges = []
    added_edges = 0
    for a_raw, b_raw in edges_raw:
        i = resolve_index(a_raw)
        j = resolve_index(b_raw)
        if i is not None and j is not None:
            fig.add_trace(go.Scatter3d(
                x=[xs[i], xs[j]],
                y=[ys[i], ys[j]],
                z=[zs[i], zs[j]],
                mode="lines",
                line=dict(width=4),
                name=f"{stars[i]['name']}–{stars[j]['name']}",
                hoverinfo="skip"
            ))
            added_edges += 1
        else:
            missing_edges.append((a_raw, b_raw))
    print(f"[INFO] edges drawn: {added_edges}")
    if missing_edges:
        print("[WARN] name mismatch in edges (not found in CSV 'name'):")
        for a, b in missing_edges:
            print(f"  - '{a}'  vs  '{b}'")

    # 軸とレイアウト
    fig.update_layout(
        title=f"{title} — 3D (Sun at origin)",
        scene=dict(
            xaxis_title="x (ly)",
            yaxis_title="y (ly)",
            zaxis_title="z (ly)",
            aspectmode="data"  # 等尺
        ),
        margin=dict(l=0, r=0, t=40, b=0)
    )

    fig.write_html(out_html, include_plotlyjs="cdn")
    print(f"[OK] Saved: {out_html}")

def parse_args(argv: List[str]):
    # 簡易オプションパーサ（--no-labels, --title）
    show_labels = True
    title = None
    files = []
    for a in argv[1:]:
        if a == "--no-labels":
            show_labels = False
        elif a.startswith("--title="):
            title = a.split("=", 1)[1]
        else:
            files.append(a)

    if len(files) not in (2, 3):
        print("Usage:")
        print("  python plot_constellation_plotly.py <stars.csv> <edges.json> <out.html> [--no-labels] [--title=...]")
        print("  python plot_constellation_plotly.py <stars.csv> <out.html> [--no-labels] [--title=...]")
        sys.exit(1)

    if len(files) == 3:
        csv_path, edges_path, out_html = files
    else:
        csv_path, out_html = files
        edges_path = None

    return csv_path, edges_path, out_html, title, show_labels

if __name__ == "__main__":
    csv_path, edges_path, out_html, title, show_labels = parse_args(sys.argv)
    plot_constellation_html(csv_path, out_html, edges_path, title=title, show_labels=show_labels)
