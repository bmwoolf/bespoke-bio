import numpy as np

from .contracts import Transform


def apply(xyz, transform: Transform):
    return (
        np.einsum("...j,ij->...i", np.asarray(xyz), np.asarray(transform.rotation))
        + transform.translation
    )


def inverse(transform: Transform):
    r = np.array(transform.rotation).T
    return Transform(rotation=r.tolist(), translation=(-r @ transform.translation).tolist())


def fit(source, destination):
    x, y = np.asarray(source), np.asarray(destination)
    if x.shape != y.shape or x.ndim != 2 or x.shape[1] != 3 or len(x) < 3:
        raise ValueError("Need at least three paired 3D points")
    if min(np.linalg.matrix_rank(x - x.mean(0)), np.linalg.matrix_rank(y - y.mean(0))) < 2:
        raise ValueError("Underdetermined alignment: collinear anchors")
    u, _, vt = np.linalg.svd((x - x.mean(0)).T @ (y - y.mean(0)))
    d = np.eye(3)
    d[-1, -1] = np.linalg.det(vt.T @ u.T)
    r = vt.T @ d @ u.T
    t = y.mean(0) - r @ x.mean(0)
    tr = Transform(rotation=r.tolist(), translation=t.tolist())
    return tr, float(np.sqrt(np.mean(np.sum((apply(x, tr) - y) ** 2, axis=1))))


def frame(origin, x_anchor, plane_anchor):
    o, a, b = map(np.asarray, (origin, x_anchor, plane_anchor))
    x = a - o
    z = np.cross(x, b - o)
    if np.linalg.norm(x) < 1e-8 or np.linalg.norm(z) < 1e-8:
        return None
    x = x / np.linalg.norm(x)
    z = z / np.linalg.norm(z)
    y = np.cross(z, x)
    return Transform(rotation=np.stack([x, y, z], axis=1).tolist(), translation=o.tolist())


def distances(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if not len(a) or not len(b):
        raise ValueError("Empty atom selection")
    return np.sqrt(np.sum((a[:, None, :] - b[None, :, :]) ** 2, axis=-1))


def pair_metrics(a, b):
    # Broad geometric triage only; not chemistry-aware energy or binding prediction.
    d = distances(a, b)
    return {
        "min_heavy_atom_distance": float(d.min()),
        "severe_clash_pairs": int((d < 1.5).sum()),
        "contact_atom_pairs": int((d < 5.0).sum()),
    }
