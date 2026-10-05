"""Fetch named arrays from the public SDRF's original ZIPs using HTTP ranges."""

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlparse

root = Path(__file__).resolve().parents[2]
local_dependencies = root / ".codex_tmp_extract" / "python_packages"
if local_dependencies.exists():
    sys.path.insert(0, str(local_dependencies))
from remotezip import RemoteZip


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("specimens", nargs="+")
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    default_cache = (Path(os.environ["LOCALAPPDATA"]) / "Codex" / "PDAC-public-data-cache"
                     if os.name == "nt" else Path.home() / ".cache" / "pdac-public-data")
    allowed = (Path(os.environ.get("PDAC_DATA_CACHE", str(default_cache))) / "emtab1791_arrays").resolve()
    if output != allowed:
        raise ValueError("Output must be the dedicated array preparation directory")
    output.mkdir(parents=True, exist_ok=True)
    with (root / "data/publication/E-MTAB-1791.sdrf.txt").open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        columns = [header.index(name) for name in
                   ("Source Name", "Derived Array Data File", "Comment [Derived ArrayExpress FTP file]")]
        records = {row[columns[0]]: (row[columns[1]], row[columns[2]]) for row in reader}

    def fetch(specimen):
        filename, archive = records[specimen]
        if Path(filename).name != filename:
            raise ValueError("Unexpected array filename")
        url = archive.replace("ftp://", "https://", 1)
        parsed = urlparse(url)
        if parsed.hostname != "ftp.ebi.ac.uk" or not parsed.path.startswith(
                "/pub/databases/microarray/data/experiment/MTAB/E-MTAB-1791/"):
            raise ValueError("Unexpected archive origin")
        target = output / filename
        partial = target.with_name(target.name + ".part")
        verified = target.with_name(target.name + ".verified.json")
        verified_partial = verified.with_name(verified.name + ".part")
        verified.unlink(missing_ok=True)
        with RemoteZip(url, timeout=(30, 180)) as archive_file:
            member = archive_file.getinfo(filename)
            with archive_file.open(member) as source, partial.open("wb") as destination:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
            if partial.stat().st_size != member.file_size:
                raise ValueError(f"Incomplete array: {specimen}")
        with partial.open("rb") as stream:
            digest = hashlib.file_digest(stream, "md5").hexdigest()
        record = {"File": filename, "Bytes": member.file_size, "MD5": digest,
                  "Archive_URL": url, "ZIP_CRC32": f"{member.CRC:08x}",
                  "Retrieval": "SDRF ZIP member (CRC checked)"}
        os.replace(partial, target)
        verified_partial.write_text(json.dumps(record), encoding="utf-8")
        os.replace(verified_partial, verified)
        print(f"Fetched {specimen} from original ZIP; CRC and size checked", flush=True)

    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(fetch, args.specimens))


if __name__ == "__main__":
    main()
