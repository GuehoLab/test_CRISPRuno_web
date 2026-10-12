#!/usr/bin/env python3
"""Rewrite the limits of the (single, defined) linear memory of a wasm module.

Usage: set_memory.py in.wasm out.wasm [--initial-pages N] [--max-pages N]
Without options it only prints the current limits. 65536 pages = 4 GiB (wasm32 limit).
Used because zig 0.13 `cc` cannot parse --max-memory values >= 1 GiB.
"""
import sys

def leb_u(b, p):
    r = s = 0
    while True:
        c = b[p]; p += 1
        r |= (c & 0x7f) << s; s += 7
        if not c & 0x80:
            return r, p

def enc_u(x):
    out = bytearray()
    while True:
        c = x & 0x7f; x >>= 7
        if x: out.append(c | 0x80)
        else:
            out.append(c); return bytes(out)

def main():
    a = sys.argv[1:]
    fin, fout = a[0], a[1] if len(a) > 1 and not a[1].startswith('--') else None
    init = int(a[a.index('--initial-pages') + 1]) if '--initial-pages' in a else None
    mx = int(a[a.index('--max-pages') + 1]) if '--max-pages' in a else None
    b = open(fin, 'rb').read()
    assert b[:4] == b'\0asm'
    p, out = 8, bytearray(b[:8])
    while p < len(b):
        sid = b[p]; size, q = leb_u(b, p + 1); body = b[q:q + size]
        if sid == 5:  # memory section
            n, r = leb_u(body, 0); assert n == 1, "expected exactly one memory"
            flags = body[r]; r += 1
            mn, r = leb_u(body, r)
            cur_max = None
            if flags & 1: cur_max, r = leb_u(body, r)
            print(f"memory: flags={flags} initial={mn} pages ({mn*65536/2**20:.0f} MiB) max={cur_max}", file=sys.stderr)
            mn = init if init is not None else mn
            new_max = mx if mx is not None else cur_max
            nb = enc_u(1) + bytes([(flags & ~1) | (1 if new_max is not None else 0)]) + enc_u(mn)
            if new_max is not None: nb += enc_u(new_max)
            body = nb + body[r:]
        out += bytes([sid]) + enc_u(len(body)) + body
        p = q + size
    if fout:
        open(fout, 'wb').write(out)

if __name__ == '__main__':
    main()
