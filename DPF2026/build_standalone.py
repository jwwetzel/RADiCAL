#!/usr/bin/env python3
"""Inline every assets/*.png into index.html as base64 data URIs -> standalone.html.
The result is ONE file that works offline on any machine (email it, USB it, open it) with
no folder dependency. Re-run after editing index.html.  Usage:  python3 build_standalone.py
"""
import base64, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "index.html"), encoding="utf-8").read()

def inline(m):
    path = m.group(1)
    full = os.path.join(HERE, path)
    if not os.path.exists(full):
        print("  WARN missing:", path); return m.group(0)
    b64 = base64.b64encode(open(full, "rb").read()).decode()
    ext = os.path.splitext(path)[1].lstrip(".").lower()
    mime = "image/png" if ext == "png" else ("image/jpeg" if ext in ("jpg","jpeg") else "image/"+ext)
    return 'src="data:%s;base64,%s"' % (mime, b64)

out = re.sub(r'src="(assets/[^"]+)"', inline, src)
out = out.replace("<title>RADiCAL — DPF 2026 (J. Wetzel)</title>",
                  "<title>RADiCAL — DPF 2026 (J. Wetzel) — standalone</title>")
dst = os.path.join(HERE, "standalone.html")
open(dst, "w", encoding="utf-8").write(out)
kb = os.path.getsize(dst) // 1024
print("wrote standalone.html  (%d KB, self-contained, %d images inlined)"
      % (kb, len(re.findall(r'data:image', out))))
