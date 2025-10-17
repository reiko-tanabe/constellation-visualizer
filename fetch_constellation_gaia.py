# fetch_constellation_gaia.py
from __future__ import annotations

import sys, csv, re, unicodedata, os, time
from typing import Any, Dict, List, Optional

from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u

def norm(s: str) -> str:
    return unicodedata.normalize("NFKC", s or "").strip()

GREEK_TO_3 = {
    "α":"alf","β":"bet","γ":"gam","δ":"del","ε":"eps","ζ":"zet","η":"eta","θ":"the",
    "ι":"iot","κ":"kap","λ":"lam","μ":"mu","ν":"nu","ξ":"xi","ο":"omi","π":"pi",
    "ρ":"rho","σ":"sig","τ":"tau","υ":"ups","φ":"phi","χ":"chi","ψ":"psi","ω":"ome",
}

def derive_simbad_id_from_name(name: str) -> Optional[str]:
    m = re.search(r"\((.+?)\)\s*$", name)
    if not m:
        return None
    token = norm(m.group(1))   # 'α Ori'
    parts = token.split()
    if len(parts) != 2:
        return None
    greek, cons = norm(parts[0]), norm(parts[1])
    code = GREEK_TO_3.get(greek)
    if not code or len(cons) < 3:
        return None
    return f"* {code} {cons}"

def ly_from_parallax_mas(parallax_mas: Optional[float]) -> Optional[float]:
    if parallax_mas is None or parallax_mas <= 0:
        return None
    d_pc = 1000.0 / parallax_mas
    return d_pc * 3.26156

def _get_num(tbl, col: str) -> Optional[float]:
    if tbl is None or col not in getattr(tbl, "colnames", []):
        return None
    v = tbl[col][0]
    try:
        if getattr(v, "mask", False):
            return None
    except Exception:
        pass
    if v is None:
        return None
    try:
        return float(v)
    except Exception:
        return None

def load_targets(csv_path: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        r = csv.DictReader(f)
        fields = [c.strip() for c in (r.fieldnames or [])]
        if "name" not in fields:
            raise ValueError(f"'name' 列が見つかりません。実フィールド: {fields}")
        for row in r:
            name = norm(row.get("name", ""))
            if not name:
                continue
            simbad_id = norm(row.get("simbad_id", "")) or derive_simbad_id_from_name(name)
            rows.append({"name": name, "simbad_id": simbad_id})
    return rows

def simbad_query_one(s: Any, ident: str, fallback_name: str) -> Dict[str, Any]:
    qid = ident or fallback_name
    tbl = s.query_object(qid)
    if tbl is None or len(tbl) == 0:
        return {
            "ra_deg": None, "dec_deg": None, "mag": None, "spec_type": None,
            "plx_mas": None, "rvz_kms": None, "raw": None
        }
    ra_deg  = _get_num(tbl, "ra")
    dec_deg = _get_num(tbl, "dec")
    Vmag    = _get_num(tbl, "V")
    sptype  = None
    if "sp_type" in tbl.colnames:
        try:
            sptype = str(tbl["sp_type"][0])
        except Exception:
            sptype = None
    # ここがポイント: plx → plx_value に名称変更
    plx_mas = _get_num(tbl, "plx_value") or _get_num(tbl, "plx")
    rvz_kms = _get_num(tbl, "rvz_radvel")
    return {
        "ra_deg": ra_deg, "dec_deg": dec_deg, "mag": Vmag, "spec_type": sptype,
        "plx_mas": plx_mas, "rvz_kms": rvz_kms, "raw": tbl
    }

def gaia_query_near_icrs(ra_deg: float, dec_deg: float) -> Optional[Dict[str, Optional[float]]]:
    """
    SIMBAD座標の近傍を Gaia DR3 で照会。
    5″で試し、ダメなら10″に拡大。各半径ごとに3回リトライ（指数バックオフ）.
    戻り値: {"parallax": float|None, "radial_velocity": float|None} or None
    """
    c = SkyCoord(ra_deg*u.deg, dec_deg*u.deg, frame="icrs")
    for arcsec in (5.0, 10.0):
        radius_deg: float = (arcsec * u.arcsec).to(u.deg).value
        query = f"""
        SELECT TOP 1
          source_id, ra, dec, parallax, radial_velocity, phot_g_mean_mag,
          distance(POINT('ICRS', ra, dec), POINT('ICRS', {c.ra.deg}, {c.dec.deg})) AS angdist
        FROM gaiadr3.gaia_source
        WHERE 1 = CONTAINS(
          POINT('ICRS', ra, dec),
          CIRCLE('ICRS', {c.ra.deg}, {c.dec.deg}, {radius_deg})
        )
        ORDER BY CASE WHEN parallax IS NULL THEN 1 ELSE 0 END ASC,
                 angdist ASC
        """
        for attempt in range(3):
            try:
                job = Gaia.launch_job_async(query)
                g = job.get_results()
                if g and len(g) > 0:
                    return {
                        "parallax": _get_num(g, "parallax"),
                        "radial_velocity": _get_num(g, "radial_velocity")
                    }
                else:
                    return {"parallax": None, "radial_velocity": None}
            except Exception as e:
                # 500などのときリトライ
                wait = 1.5 * (2 ** attempt)  # 1.5s, 3s, 6s
                print(f"[WARN] Gaia failed (try {attempt+1}/3, {arcsec}\" ): {e} -> retry in {wait:.1f}s")
                time.sleep(wait)
        # 半径を広げて再挑戦
    return None

def main():
    in_csv  = sys.argv[1] if len(sys.argv) > 1 else "Orion.csv"
    out_csv = sys.argv[2] if len(sys.argv) > 2 else "Orion_gaia.csv"

    if not os.path.exists(in_csv):
        raise FileNotFoundError(
            f"Input CSV not found: {in_csv}\n"
            f"Current dir: {os.getcwd()}\n"
            f"Hint: pass a relative path like 'orion/Orion.csv'"
        )

    targets = load_targets(in_csv)

    # SIMBAD（plx→plx_value に合わせる）
    s = Simbad()
    s.add_votable_fields("ra", "dec", "sp", "V", "plx_value", "rvz_radvel")

    out_rows: List[Dict[str, Any]] = []

    for t in targets:
        name = t["name"]
        simbad_id = t["simbad_id"]
        print(f"\n[FETCH] {name}   (SIMBAD: {simbad_id or name})")

        row: Dict[str, Any] = {
            "name": name, "simbad_id": simbad_id,
            "ra_deg": None, "dec_deg": None, "dist_ly": None,
            "mag": None, "spec_type": None,
            "parallax_mas": None, "radial_velocity": None
        }

        # --- SIMBAD
        try:
            r = simbad_query_one(s, simbad_id, name)
            row["ra_deg"] = r["ra_deg"]
            row["dec_deg"] = r["dec_deg"]
            row["mag"] = r["mag"]
            row["spec_type"] = r["spec_type"]
            simbad_plx = r["plx_mas"]
            simbad_rvz = r["rvz_kms"]
        except Exception as e:
            print(f"[WARN] SIMBAD failed: {e}")
            simbad_plx = None
            simbad_rvz = None

        # --- Gaia（近傍照会 + リトライ＆半径拡大）
        if row["ra_deg"] is not None and row["dec_deg"] is not None:
            g = gaia_query_near_icrs(row["ra_deg"], row["dec_deg"])
            if g is not None:
                row["parallax_mas"] = g["parallax"]
                row["radial_velocity"] = g["radial_velocity"]
                row["dist_ly"] = ly_from_parallax_mas(g["parallax"])

        # --- 補完：Gaiaで取れなかったらSIMBADの値で後詰め
        if row["parallax_mas"] is None and simbad_plx is not None:
            row["parallax_mas"] = simbad_plx
            row["dist_ly"] = ly_from_parallax_mas(simbad_plx)
        if row["radial_velocity"] is None and simbad_rvz is not None:
            row["radial_velocity"] = simbad_rvz

        out_rows.append(row)

    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "name","simbad_id","ra_deg","dec_deg","dist_ly","mag","spec_type","parallax_mas","radial_velocity"
            ]
        )
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    print(f"\n[OK] saved: {out_csv}")

if __name__ == "__main__":
    main()
