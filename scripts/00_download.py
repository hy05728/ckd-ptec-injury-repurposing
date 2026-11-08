#!/usr/bin/env python3
"""00_download.py — 用 urllib 稳健下载 checkpoint 与 GSE131882 原始数据（带断点续传）。"""
import os, sys, time, urllib.request

FILES = [
    ("https://zenodo.org/api/records/17437512/files/model_fold_0.pt/content",
     r"E:/workbuddy/2026-11-08-00-03-12/data/checkpoints/model_fold_0.pt"),
    ("https://zenodo.org/api/records/17437512/files/model_fold_1.pt/content",
     r"E:/workbuddy/2026-11-08-00-03-12/data/checkpoints/model_fold_1.pt"),
    ("https://zenodo.org/api/records/17437512/files/model_fold_2.pt/content",
     r"E:/workbuddy/2026-11-08-00-03-12/data/checkpoints/model_fold_2.pt"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131882/suppl/GSE131882_RAW.tar",
     r"E:/workbuddy/2026-11-08-00-03-12/data/GSE131882_RAW.tar"),
]

def download(url, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    existing = os.path.getsize(dest) if os.path.exists(dest) else 0
    headers = {"Range": f"bytes={existing}-"} if existing else {}
    req = urllib.request.Request(url, headers=headers)
    mode = "ab" if existing else "wb"
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, mode) as f:
        clen = r.headers.get("Content-Length")
        total = int(clen) if clen else None
        n = 0
        while True:
            chunk = r.read(1024 * 256)
            if not chunk:
                break
            f.write(chunk)
            n += len(chunk)
            got = existing + n
            if total:
                pct = got / (existing + total) * 100
                print(f"  {os.path.basename(dest)}: {got/(1<<20):.1f} MB  ({pct:.0f}%)", end="\r")
            else:
                print(f"  {os.path.basename(dest)}: {got/(1<<20):.1f} MB", end="\r")
    print(f"\nDONE {os.path.basename(dest)} -> {os.path.getsize(dest)} bytes")

for url, dest in FILES:
    print(f"[dl] {os.path.basename(dest)}")
    for attempt in range(3):
        try:
            download(url, dest)
            break
        except Exception as e:
            print(f"  attempt {attempt+1} failed: {type(e).__name__}: {e}")
            time.sleep(8)
    else:
        print(f"FAILED {dest}")
print("ALL DONE")
