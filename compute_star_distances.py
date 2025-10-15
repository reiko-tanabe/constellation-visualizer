import csv, math, os, sys, json, unicodedata, re
from typing import Any, Dict, List, Optional, Tuple

# ---------- helpers (同じロジックでCSVを読む/名前を揃える) ----------
def norm(s: str) -> str:
    if s is None: return ""
    return unicodedata.normalize("NFKC", s).strip()

def base_key(s: str) -> str:
    s2 = norm(s)
    s2 = re.sub(r"\s*\(.*?\)\s*$", "", s2)  # 末尾の括弧表記を削除
    s2 = re.sub(r"\s+", " ", s2)
    return s2.lower()

def _float_or_none(v: Any) -> Optional[float]:
    try:
        if v is None: return None
        s = str(v).strip()
        if s == "" or s.lower() in ("nan", "none"): return None
        return float(s)
    except Exception:
        return None

def sph_to_cart(ra_deg: float, dec_deg: float, dist_ly: float) -> Tuple[float,float,float]:
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    x = dist_ly * math.cos(dec) * math.cos(ra)
    y = dist_ly * math.cos(dec) * math.sin(ra)
    z = dist_ly * math.sin(dec)
    return x, y, z

def load_stars(csv_path: str) -> List[Dict[str, Any]]:
    """ra_deg/dec_deg/dist_ly が入ったCSV（utf-8 / utf-8-sig）を読み、x,y,zを付与して返す"""
    stars: List[Dict[str, Any]] = []
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        r = csv.DictReader(f)
        fields = [c.strip() for c in (r.fieldnames or [])]
        need = {"name","ra_deg","dec_deg","dist_ly"}
        if not need.issubset(set(fields)):
            raise ValueError(f"CSVに {need} が必要です。実ヘッダ: {fields}")

        for row in r:
            name = row.get("name","").strip()
            ra   = _float_or_none(row.get("ra_deg"))
            dec  = _float_or_none(row.get("dec_deg"))
            dist = _float_or_none(row.get("dist_ly"))
            if name == "" or ra is None or dec is None or dist is None or dist <= 0:
                # 計算できない行はスキップ
                continue
            x,y,z = sph_to_cart(ra, dec, dist)
            stars.append({
                "name": name,
                "name_norm": norm(name),
                "name_base": base_key(name),
                "ra_deg": ra, "dec_deg": dec, "dist_ly": dist,
                "x": x, "y": y, "z": z
            })
    return stars

def load_edges(edges_path: Optional[str]) -> List[Tuple[str,str]]:
    if not edges_path: return []
    if not os.path.exists(edges_path): return []
    with open(edges_path, encoding="utf-8-sig") as f:
        raw = json.load(f)
    out: List[Tuple[str,str]] = []
    for item in raw:
        if isinstance(item, (list, tuple)) and len(item) == 2:
            out.append((str(item[0]), str(item[1])))
    return out

def resolve_index(name_raw: str, name_to_idx: Dict[str,int], base_to_idx: Dict[str,int]) -> Optional[int]:
    n = norm(name_raw)
    if n in name_to_idx: return name_to_idx[n]
    b = base_key(name_raw)
    if b in base_to_idx: return base_to_idx[b]
    return None

def dist3(a, b) -> float:
    dx = a["x"] - b["x"]
    dy = a["y"] - b["y"]
    dz = a["z"] - b["z"]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

# ---------- main ----------
def main(argv: List[str]) -> None:
    if len(argv) < 2:
        print("Usage:")
        print("  python compute_star_distances.py <stars.csv> [edges.json] [--all-pairs] [--out=distances.csv]")
        sys.exit(1)

    csv_path = argv[1]
    edges_path = None
    out_csv = "distances.csv"
    all_pairs = False

    for a in argv[2:]:
        if a.endswith(".json"):
            edges_path = a
        elif a == "--all-pairs":
            all_pairs = True
        elif a.startswith("--out="):
            out_csv = a.split("=",1)[1]

    stars = load_stars(csv_path)
    if not stars:
        print("[ERROR] no stars with valid (ra, dec, dist_ly)")
        sys.exit(2)

    # 名前→index辞書（完全一致＆括弧落としキー）
    name_to_idx = {s["name_norm"]: i for i,s in enumerate(stars)}
    # baseキーは重複があれば曖昧なので除外
    tmp = {}
    base_to_idx = {}
    for i,s in enumerate(stars):
        b = s["name_base"]
        tmp[b] = tmp.get(b,0)+1
        if b not in base_to_idx:
            base_to_idx[b] = i
    for b,cnt in tmp.items():
        if cnt>1 and b in base_to_idx:
            del base_to_idx[b]

    rows_out: List[Dict[str,Any]] = []

    if all_pairs:
        # 全組み合わせ
        for i in range(len(stars)):
            for j in range(i+1, len(stars)):
                a, b = stars[i], stars[j]
                d = dist3(a,b)
                rows_out.append({
                    "star_a": a["name"],
                    "star_b": b["name"],
                    "distance_ly": d
                })
        print(f"[INFO] pairs computed: {len(rows_out)}")

    else:
        # edges.json が必要
        edges = load_edges(edges_path)
        print(f"[INFO] edges loaded: {len(edges)}")
        missing = 0
        for a_raw, b_raw in edges:
            ia = resolve_index(a_raw, name_to_idx, base_to_idx)
            ib = resolve_index(b_raw, name_to_idx, base_to_idx)
            if ia is None or ib is None:
                print(f"[WARN] name mismatch: '{a_raw}' vs '{b_raw}'")
                missing += 1
                continue
            a, b = stars[ia], stars[ib]
            rows_out.append({
                "star_a": a["name"],
                "star_b": b["name"],
                "distance_ly": dist3(a,b)
            })
        print(f"[INFO] edges resolved: {len(rows_out)}  (missing: {missing})")

    # 出力
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["star_a","star_b","distance_ly"])
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"[OK] saved: {out_csv}")

if __name__ == "__main__":
    main(sys.argv)
