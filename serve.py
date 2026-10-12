#!/usr/bin/env python3
"""Serve the web/ folder on http://localhost:8000 (OPFS needs a secure context; localhost counts).
Sends COOP/COEP headers so that the pages are cross-origin isolated (needed for multi-threading).

    python serve.py [port]
"""
import http.server, os, sys, functools

class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map,
                      '.wasm': 'application/wasm', '.js': 'text/javascript'}
    def end_headers(self):
        self.send_header('Cache-Control', 'no-cache')
        # cross-origin isolation -> SharedArrayBuffer -> multi-threaded build (minibwa64-mt.js)
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        super().end_headers()

port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'web')
print(f'Serving {root} at http://localhost:{port}/')
http.server.ThreadingHTTPServer(('127.0.0.1', port), functools.partial(Handler, directory=root)).serve_forever()
