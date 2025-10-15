# pick_landmarks.py — Mac / QtAgg 用 安定版
import argparse, json, csv
import matplotlib
# バックエンドを QtAgg に固定（macOSX バックエンドだとエラーが出やすい）
matplotlib.use("QtAgg", force=True)
import matplotlib.pyplot as plt
import numpy as np

def read_star_names(csv_path):
    """CSVから星の名前を読み取る"""
    names = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            names.append(row["name"])
    return names

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("image_png")
    ap.add_argument("stars_csv")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    names = read_star_names(args.stars_csv)
    img = plt.imread(args.image_png)
    H, W = img.shape[:2]

    out_json = args.out or (args.image_png.rsplit(".",1)[0] + "_landmarks.json")
    print(f"[INFO] will save to: {out_json}")

    # 既存があれば読み込み
    points = {}
    try:
        with open(out_json, "r", encoding="utf-8") as f:
            points = json.load(f)
            print(f"[INFO] loaded existing: {out_json}")
    except Exception:
        pass

    # 未入力リスト
    order = [n for n in names if n not in points]

    fig, ax = plt.subplots(figsize=(min(10, W/100), min(10, H/100)))
    fig.canvas.manager.set_window_title("pick_landmarks — Left click=record, u=undo, s=save, q=quit")
    ax.imshow(img)
    dots = ax.scatter([], [], s=40, c="cyan")
    tip  = ax.text(5, 10, "", color="yellow", fontsize=12, va="top")
    ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout()

    history = []  # (name, x, y)

    def refresh():
        xs, ys = [], []
        for n in names:
            if n in points:
                x, y = points[n]; xs.append(x); ys.append(y)
        offset = np.column_stack([xs, ys]) if xs else np.empty((0,2))
        dots.set_offsets(offset)
        next_name = order[0] if order else "(done)"
        tip.set_text(f"Next: {next_name}   saved: {len(points)}/{len(names)}")
        fig.canvas.draw_idle()

    def onclick(ev):
        if ev.inaxes != ax:
            return
        if not order:
            return

        # データ座標を直接取得（これが安定）
        if ev.xdata is None or ev.ydata is None:
            return
        x_data, y_data = float(ev.xdata), float(ev.ydata)

        name = order.pop(0)
        points[name] = [x_data, y_data]
        history.append((name, x_data, y_data))
        print(f"[CLICK] {name}: data=({x_data:.1f}, {y_data:.1f})")

        # 自動保存
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(points, f, ensure_ascii=False, indent=2)
        print(f"[AUTO-SAVE] {out_json}")
        refresh()

    def onkey(ev):
        if ev.key in ("u", "backspace"):
            if history:
                name, _, _ = history.pop()
                points.pop(name, None)
                order.insert(0, name)
                print(f"[UNDO] {name}")
                refresh()
        elif ev.key == "s":
            with open(out_json, "w", encoding="utf-8") as f:
                json.dump(points, f, ensure_ascii=False, indent=2)
            print(f"[SAVE] {out_json}")
        elif ev.key == "q":
            plt.close(fig)

    def onmove(ev):
        if ev.inaxes != ax:
            return
        if ev.xdata is None or ev.ydata is None:
            return
        x_data, y_data = float(ev.xdata), float(ev.ydata)
        tip.set_text(f"Next: {order[0] if order else '(done)'}   "
                     f"saved: {len(points)}/{len(names)}   "
                     f"x={x_data:.1f}, y={y_data:.1f}")
        fig.canvas.draw_idle()

    # イベント登録
    fig.canvas.mpl_connect("button_press_event", onclick)
    fig.canvas.mpl_connect("key_press_event", onkey)
    fig.canvas.mpl_connect("motion_notify_event", onmove)

    refresh()
    plt.show(block=True)   # ← ウィンドウを閉じるまで待機

    # 閉じる際に自動保存
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(points, f, ensure_ascii=False, indent=2)
    print(f"[EXIT SAVE] {out_json}")

if __name__ == "__main__":
    main()
