"""
fetch_aries_gaia.py — 牡羊座主要恒星の位置と距離をSIMBAD・Gaia・HIP2から取得（最小・完全版）
目的: 3D星座描画に必要な最小情報 (RA, Dec, 距離, 等級, 分類)
"""

import re, time
import pandas as pd
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord

# ============ 手動補完（最終手段） ============
MANUAL_OVERRIDES = {
    "Mesarthim (γ Ari)": {"dist_ly": 164.0},  # 文献値に基づく目安
}

# ============ 主要星リスト ============
STAR_CANDIDATES = {
    "Hamal (α Ari)": "* alf Ari",
    "Sheratan (β Ari)": "* bet Ari",
    "Mesarthim (γ Ari)": "* gam1 Ari",
    "Eta Arietis (η Ari)": "* eta Ari",
}

# ============ SIMBAD設定 ============
simbad = Simbad()
simbad.add_votable_fields("ra","dec","sp","V","G","pmra","pmdec","parallax","ids","rvz_radvel")

# ============ HIP2設定 ============
Vizier.ROW_LIMIT = 5
def fetch_hip2(hip):
    """HIP2 (I/311/hip2) から視差・PM取得"""
    if hip is None:
        return None
    try:
        res = Vizier(columns=["HIP","Plx","pmRA","pmDE"]).query_constraints(
            catalog="I/311/hip2", HIP=str(hip)
        )
        if len(res) and len(res[0]):
            r = res[0][0]
            return {
                "parallax_mas": float(r["Plx"]) if r["Plx"] else None,
                "pmra": float(r["pmRA"]) if r["pmRA"] else None,
                "pmdec": float(r["pmDE"]) if r["pmDE"] else None,
            }
    except Exception:
        pass
    return None

def get_hip_from_simbad_objectids(simbad_id: str):
    """objectids から HIP 番号取得"""
    try:
        ids_tbl = Simbad().query_objectids(simbad_id)
        if ids_tbl is None:
            return None
        for row in ids_tbl:
            oid = str(row["id"])
            if oid.startswith("HIP "):
                return int(oid.split()[1])
    except Exception:
        pass
    return None

# ============ ID解析 ============
def parse_ids_to_gaia_hip(ids_str: str):
    gaia_id = hip = None
    if ids_str:
        m = re.search(r"Gaia DR3\s+(\d+)", ids_str)
        if m: gaia_id = int(m.group(1))
        m2 = re.search(r"HIP\s+(\d+)", ids_str)
        if m2: hip = int(m2.group(1))
    return gaia_id, hip

# ============ Gaia検索 ============
def fetch_gaia_by_source_id(source_id):
    q = f"""
      SELECT source_id, ra, dec, parallax, phot_g_mean_mag, radial_velocity
      FROM gaiadr3.gaia_source WHERE source_id = {source_id}
    """
    job = Gaia.launch_job(q)
    r = job.get_results()
    return r[0] if len(r) else None

def fetch_gaia_by_hip(hip):
    q = f"""
      SELECT g.source_id, g.ra, g.dec, g.parallax, g.phot_g_mean_mag, g.radial_velocity
      FROM gaiadr3.hipparcos2_best_neighbour AS h
      JOIN gaiadr3.gaia_source AS g ON g.source_id = h.source_id
      WHERE h.original_ext_source_id = {hip}
    """
    job = Gaia.launch_job(q)
    r = job.get_results()
    return r[0] if len(r) else None

def fetch_gaia_by_cone(ra_deg, dec_deg):
    for radius_arcsec in (5.0, 15.0, 30.0):
        q = f"""
          SELECT TOP 1 source_id, ra, dec, parallax, phot_g_mean_mag, radial_velocity
          FROM gaiadr3.gaia_source
          WHERE 1=CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra_deg}, {dec_deg}, {radius_arcsec}/3600.0))
        """
        try:
            job = Gaia.launch_job(q)
            r = job.get_results()
            if len(r):
                r[0]["_search_radius_arcsec"] = radius_arcsec
                return r[0]
        except Exception:
            time.sleep(1)
            continue
    return None

def gaia_pick(row, col):
    try:
        v = row[col]
        if hasattr(v, "mask"):
            return float(v) if not v.mask else None
        return float(v)
    except Exception:
        return None

# ============ メイン ============
rows = []
for disp_name, simbad_id in STAR_CANDIDATES.items():
    print(f"\n[FETCH] {disp_name}")
    tbl = simbad.query_object(simbad_id)
    if tbl is None or len(tbl) == 0:
        print(f"[WARN] SIMBAD not found: {disp_name}")
        continue

    ra = float(tbl["ra"][0])
    dec = float(tbl["dec"][0])
    vmag = float(tbl["V"][0]) if "V" in tbl.colnames and tbl["V"][0] is not None else None
    g_from_simbad = float(tbl["G"][0]) if "G" in tbl.colnames and tbl["G"][0] is not None else None
    sp = tbl["sp_type"][0] if "sp_type" in tbl.colnames else None
    plx_simbad = float(tbl["parallax"][0]) if "parallax" in tbl.colnames and tbl["parallax"][0] is not None else None
    ids_str = str(tbl["ids"][0]) if "ids" in tbl.colnames and tbl["ids"][0] is not None else ""
    rv_simbad = float(tbl["rvz_radvel"][0]) if "rvz_radvel" in tbl.colnames and tbl["rvz_radvel"][0] is not None else None

    gaia_id, hip = parse_ids_to_gaia_hip(ids_str)
    if hip is None:
        hip = get_hip_from_simbad_objectids(simbad_id)

    g = None
    if gaia_id: g = fetch_gaia_by_source_id(gaia_id)
    if g is None and hip: g = fetch_gaia_by_hip(hip)
    if g is None: g = fetch_gaia_by_cone(ra, dec)

    parallax_mas = dist_ly = gmag = rv = None
    if g is not None:
        parallax_mas = gaia_pick(g, "parallax")
        gmag = gaia_pick(g, "phot_g_mean_mag")
        rv = gaia_pick(g, "radial_velocity")

    # HIP2補完
    hip2 = fetch_hip2(hip)
    if hip2 and parallax_mas is None and hip2["parallax_mas"]:
        parallax_mas = hip2["parallax_mas"]

    # 視差→距離
    parallax_chain = [parallax_mas, plx_simbad, hip2["parallax_mas"] if hip2 else None]
    plx_used = next((p for p in parallax_chain if p is not None and p > 0), None)
    if plx_used:
        dist_ly = (1000.0 / plx_used) * 3.26156

    # RVフォールバック
    if rv is None and rv_simbad is not None:
        rv = rv_simbad

    # 等級フォールバック
    mag_out = (
        vmag
        if vmag is not None
        else (g_from_simbad if g_from_simbad is not None else gmag)
    )

    # 手動補完
    if dist_ly is None and disp_name in MANUAL_OVERRIDES and "dist_ly" in MANUAL_OVERRIDES[disp_name]:
        dist_ly = float(MANUAL_OVERRIDES[disp_name]["dist_ly"])

    rows.append(dict(
        name=disp_name,
        ra_deg=ra, dec_deg=dec,
        dist_ly=dist_ly,
        mag=mag_out,
        spec_type=sp,
        parallax_mas=plx_used,
        radial_velocity=rv
    ))

# ============ 出力 ============
df = pd.DataFrame(rows)
out_csv = "Aries_gaia.csv"
df.to_csv(out_csv, index=False, encoding="utf-8-sig")
print(f"\n[OK] saved: {out_csv}")
print(df)
