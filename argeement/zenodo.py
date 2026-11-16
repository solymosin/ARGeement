"""Publication of a release on Zenodo, as a new version of a versioned dataset record.

    python -m argeement zenodo --package <package dir> [--config zenodo.yaml] [--sandbox] [--dry-run]
    python -m argeement zenodo --release <release dir> ...     (packages the release first)

Steps: new version of the record, removal of the previous files from the draft, upload of the
archive and its checksum, metadata (version, date, database versions), publication. The DOIs are
written to zenodo.json in the package (or release) directory; argeement site reads them from there.

The first version is created with --new (once, when zenodo.yaml has no id yet); its id must then be
written to deposition_id (sandbox: sandbox_deposition_id) in zenodo.yaml.

The access token is read from the ZENODO_TOKEN environment variable or from
~/.config/argeement/zenodo_token (sandbox: ZENODO_SANDBOX_TOKEN or zenodo_sandbox_token); it is never
stored in the repository. The id of the latest published version is kept in zenodo_state.json
(ignored by git).
"""
import datetime as dt
import json
import os
import urllib.error
import urllib.request
from pathlib import Path

import yaml

from .package import package
from .terms import DATA_LICENCE, licence_html

API = {"zenodo": "https://zenodo.org/api", "sandbox": "https://sandbox.zenodo.org/api"}
NAME = {"amrfinder": "AMRFinderPlus", "card": "CARD", "resfinder": "ResFinder",
        "disinfinder": "DisinFinder", "pointfinder": "PointFinder"}


# an explicit client name: some web application firewalls reject the default "Python-urllib/x.y"
USER_AGENT = "ARGeement (+https://github.com/solymosin/ARGeement)"


def _token(sandbox):
    env = "ZENODO_SANDBOX_TOKEN" if sandbox else "ZENODO_TOKEN"
    if os.environ.get(env):
        return os.environ[env]
    f = Path.home() / ".config" / "argeement" / ("zenodo_sandbox_token" if sandbox else "zenodo_token")
    if f.exists():
        return f.read_text().strip()
    raise SystemExit(f"No Zenodo token: set the {env} environment variable or create {f}")


class Zenodo:
    def __init__(self, token, sandbox=False):
        # ZENODO_API: for testing only (local mock server)
        self.base = os.environ.get("ZENODO_API") or API["sandbox" if sandbox else "zenodo"]
        self.token = token

    def call(self, method, url, data=None, raw=None, ctype="application/json"):
        url = url if url.startswith("http") else self.base + url
        body = raw if raw is not None else (json.dumps(data).encode() if data is not None else None)
        req = urllib.request.Request(url, data=body, method=method)
        req.add_header("Authorization", f"Bearer {self.token}")
        req.add_header("User-Agent", USER_AGENT)
        req.add_header("Accept", "application/json")
        if body is not None:
            req.add_header("Content-Type", ctype)
        try:
            with urllib.request.urlopen(req, timeout=600) as r:
                txt = r.read().decode()
                return json.loads(txt) if txt else {}
        except urllib.error.HTTPError as e:
            raise SystemExit(f"Zenodo error ({method} {url}): {e.code} {e.read().decode()[:800]}")

    def upload(self, bucket, path):
        size = path.stat().st_size
        with open(path, "rb") as fh:
            req = urllib.request.Request(f"{bucket}/{path.name}", data=fh, method="PUT")
            req.add_header("Authorization", f"Bearer {self.token}")
            req.add_header("User-Agent", USER_AGENT)
            req.add_header("Content-Type", "application/octet-stream")
            req.add_header("Content-Length", str(size))
            try:
                with urllib.request.urlopen(req, timeout=3600) as r:
                    return json.loads(r.read().decode())
            except urllib.error.HTTPError as e:
                raise SystemExit(f"Upload error: {e.code} {e.read().decode()[:800]}")


def metadata(cfg, man, release_name):
    vers = man.get("database_versions", {})
    vtab = "".join(f"<li>{NAME.get(k, k)}: {v}</li>" for k, v in sorted(vers.items()))
    desc = (cfg.get("description", "").strip()
            + f"<p>Release <b>{release_name}</b> (profile <code>{man.get('profile', '')}</code>), "
              f"created {man.get('created', '')[:10]}. Database versions:</p><ul>{vtab}</ul>"
            + "<h3>Licence</h3>" + licence_html())
    md = dict(title=f"{cfg['title']} – {release_name}", upload_type="dataset",
              description=desc, creators=cfg["creators"], version=release_name,
              publication_date=dt.date.today().isoformat(),
              keywords=cfg.get("keywords", []), license=cfg.get("license", DATA_LICENCE.lower()))
    for k in ("related_identifiers", "communities", "notes", "access_right"):
        if cfg.get(k):
            md[k] = cfg[k]
    return md


def publish(release=None, config="zenodo.yaml", sandbox=False, dry_run=False, state_file=None,
            package_dir=None, new_record=False):
    """release: release directory (the package is made here); package_dir: a finished package
    directory, possibly downloaded to another machine (archive, .sha256, manifest.json, package.json)."""
    if package_dir:
        pkd = Path(package_dir).resolve()
        pj = json.loads((pkd / "package.json").read_text())
        man = json.loads((pkd / "manifest.json").read_text())
        rel_name, out_dir = pj["run"], pkd
        archive = pkd / pj["archive"]
        got = __import__("hashlib").sha256(archive.read_bytes()).hexdigest()
        if got != pj["sha256"]:
            raise SystemExit(f"Checksum mismatch for {archive.name}: incomplete or corrupted download?")
    else:
        rel = Path(release).resolve()
        man = json.loads((rel / "manifest.json").read_text())
        rel_name, out_dir, archive = rel.name, rel, None
    cfg = yaml.safe_load(Path(config).read_text())
    state_file = Path(state_file or Path(config).with_name(
        "zenodo_state_sandbox.json" if sandbox else "zenodo_state.json"))
    state = json.loads(state_file.read_text()) if state_file.exists() else {}
    dep_id = state.get("deposition_id") or (cfg.get("sandbox_deposition_id") if sandbox
                                            else cfg.get("deposition_id"))
    key = "sandbox_deposition_id" if sandbox else "deposition_id"
    if dep_id and new_record:
        raise SystemExit(f"--new: the record already exists ({key}: {dep_id}); every further release "
                         "becomes a new version of it, run without --new")
    if not dep_id and not new_record:
        raise SystemExit(f"No Zenodo record yet ({key} is empty in {config}). To create the first "
                         "version of a new record, run once with --new (publish_local.sh: NEW=1).")
    if state.get("published", {}).get(rel_name):
        res = state["published"][rel_name]
        print(f"Release {rel_name} is already published: {res['doi']}")
        (out_dir / "zenodo.json").write_text(json.dumps(res, indent=2))
        return res

    if archive is None:
        archive = Path(package(rel, rel / "package")["archive"])
    md = metadata(cfg, man, rel_name)
    if dry_run:
        print(json.dumps(dict(deposition_id=dep_id or "(new record)", archive=str(archive), metadata=md), indent=2,
                         ensure_ascii=False))
        return None

    z = Zenodo(_token(sandbox), sandbox)
    if dep_id:
        nv = z.call("POST", f"/deposit/depositions/{dep_id}/actions/newversion")
        draft = z.call("GET", nv["links"]["latest_draft"])
        for f in draft.get("files", []):                   # files of the previous version are not inherited
            z.call("DELETE", f["links"]["self"])
    else:                                                  # first version of a new record
        draft = z.call("POST", "/deposit/depositions", data={})
    z.upload(draft["links"]["bucket"], archive)
    sha = Path(f"{archive}.sha256")
    if sha.exists():
        z.upload(draft["links"]["bucket"], sha)
    z.call("PUT", f"/deposit/depositions/{draft['id']}", data={"metadata": md})
    pub = z.call("POST", f"/deposit/depositions/{draft['id']}/actions/publish")
    res = dict(release=rel_name, doi=pub["doi"], conceptdoi=pub.get("conceptdoi", ""),
               record_url=pub["links"].get("record_html") or pub["links"].get("html", ""),
               deposition_id=pub["id"], sandbox=sandbox,
               published=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    (out_dir / "zenodo.json").write_text(json.dumps(res, indent=2))
    state["deposition_id"] = pub["id"]
    state.setdefault("published", {})[rel_name] = res
    state_file.write_text(json.dumps(state, indent=2))
    print(f"Published: doi:{res['doi']} (concept: {res['conceptdoi']})")
    if not dep_id:
        print(f"First version of a new record. Write its id to {config}:\n  {key}: {pub['id']}")
    return res
