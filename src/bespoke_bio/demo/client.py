"""Bounded HTTP transport. Tokens remain on the server; redirects are not followed."""

import json
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener


class ServiceError(Exception):
    def __init__(self, status, message):
        self.status, self.message = status, message
        super().__init__(message)


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


class DesignService:
    def __init__(self, base_url, token=None):
        parsed = urlsplit(base_url)
        if (
            not parsed.hostname
            or parsed.username
            or parsed.password
            or parsed.query
            or parsed.fragment
            or parsed.path not in {"", "/"}
            or (
                parsed.scheme != "https"
                and not (
                    parsed.scheme == "http" and parsed.hostname in {"127.0.0.1", "localhost", "::1"}
                )
            )
        ):
            raise ValueError("Design API must use an HTTPS origin or local HTTP")
        self.base_url, self.token = base_url.rstrip("/"), token
        self.opener = build_opener(NoRedirect())

    def request(self, path, payload=None, *, binary=False):
        if not path.startswith("/v1/") or ".." in path or "?" in path or "#" in path:
            raise ValueError("Invalid design API path")
        headers = {"Accept": "application/json"}
        if self.token:
            headers["Authorization"] = "Bearer " + self.token
        data = None
        if payload is not None:
            headers["Content-Type"] = "application/json"
            data = json.dumps(payload, allow_nan=False).encode()
        try:
            with self.opener.open(
                Request(self.base_url + path, data, headers), timeout=60
            ) as response:
                limit = 128 * 1024 * 1024 if binary else 8 * 1024 * 1024
                body = response.read(limit + 1)
                if len(body) > limit:
                    raise ServiceError(502, "Design service response exceeds the download limit")
        except HTTPError as exc:
            messages = {
                401: "Design service authentication failed; configure the service token",
                403: "Design service denied this request",
                404: "Design service job or artifact was not found",
                409: "Results are not ready, or this job cannot perform the requested action",
                422: "Design service rejected the input; inspect its documented contract",
            }
            # Never forward arbitrary upstream bodies that could contain credentials or paths.
            raise ServiceError(
                exc.code if exc.code in messages else 502,
                messages.get(exc.code, "Design service request failed"),
            ) from exc
        except (URLError, TimeoutError, OSError) as exc:
            raise ServiceError(
                502, "Design service could not be reached; submission may be unknown"
            ) from exc
        if binary:
            return body
        try:
            value = json.loads(body)
            if not isinstance(value, dict):
                raise TypeError("Expected a JSON object")
            return value
        except (ValueError, TypeError, UnicodeDecodeError) as exc:
            raise ServiceError(502, "Design service returned invalid JSON") from exc
