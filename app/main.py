"""HTTP API for concatemer decoding.

Endpoints
---------
POST /api/concatemers/decode
    Body: {"reference", "read", "copies", "max_edits_per_copy"}
GET  /healthz
    Liveness/readiness probe.

Only the Python standard library is used.
"""

from __future__ import annotations

import json
import os
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Dict

from .solver import ConstraintFailure, InvalidRequest, decode

MAX_BODY_BYTES = 64 * 1024


class Handler(BaseHTTPRequestHandler):
    server_version = "ConcatemerDecoder/1.0"

    def _send_json(self, status: int, payload: Dict[str, Any]) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:  # noqa: N802 - stdlib hook
        if self.path.split("?", 1)[0] in ("/healthz", "/health"):
            self._send_json(HTTPStatus.OK, {"status": "ok"})
            return
        self._send_json(HTTPStatus.NOT_FOUND, {
            "error": "not_found",
            "message": f"Unknown path: {self.path}",
        })

    def do_POST(self) -> None:  # noqa: N802 - stdlib hook
        if self.path.split("?", 1)[0] != "/api/concatemers/decode":
            self._send_json(HTTPStatus.NOT_FOUND, {
                "error": "not_found",
                "message": f"Unknown path: {self.path}",
            })
            return

        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0:
            self._send_json(HTTPStatus.BAD_REQUEST, {
                "error": "bad_request",
                "message": "Expected a JSON request body.",
                "field": "$",
            })
            return
        if length > MAX_BODY_BYTES:
            self._send_json(HTTPStatus.REQUEST_ENTITY_TOO_LARGE, {
                "error": "bad_request",
                "message": "Request body too large.",
                "field": "$",
            })
            return

        raw = self.rfile.read(length)
        try:
            data = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            self._send_json(HTTPStatus.BAD_REQUEST, {
                "error": "bad_request",
                "message": f"Body is not valid JSON: {exc}",
                "field": "$",
            })
            return

        if not isinstance(data, dict):
            self._send_json(HTTPStatus.BAD_REQUEST, {
                "error": "bad_request",
                "message": "Request body must be a JSON object.",
                "field": "$",
            })
            return

        required = ("reference", "read", "copies", "max_edits_per_copy")
        missing = [name for name in required if name not in data]
        if missing:
            self._send_json(HTTPStatus.BAD_REQUEST, {
                "error": "bad_request",
                "message": f"Missing required field(s): {', '.join(missing)}",
                "field": missing[0],
            })
            return

        try:
            result = decode(
                data["reference"],
                data["read"],
                data["copies"],
                data["max_edits_per_copy"],
            )
        except InvalidRequest as exc:
            self._send_json(HTTPStatus.UNPROCESSABLE_ENTITY, {
                "error": "invalid_request",
                "message": exc.message,
                "field": exc.field,
            })
            return
        except ConstraintFailure as exc:
            payload = dict(exc.payload)
            payload.setdefault("endpoint", "/api/concatemers/decode")
            self._send_json(HTTPStatus.UNPROCESSABLE_ENTITY, payload)
            return

        self._send_json(HTTPStatus.OK, result)

    def log_message(self, fmt: str, *args: Any) -> None:
        if os.environ.get("QUIET"):
            return
        super().log_message(fmt, *args)


def create_server(host: str, port: int) -> ThreadingHTTPServer:
    server = ThreadingHTTPServer((host, port), Handler)
    server.daemon_threads = True
    return server


def main() -> None:
    host = os.environ.get("HOST", "0.0.0.0")
    port = int(os.environ.get("PORT", "8080"))
    server = create_server(host, port)
    print(f"concatemer decoder listening on {host}:{port}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
