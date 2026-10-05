"""Check completed-file integrity without contacting the archive server."""

import csv
import hashlib
import importlib.util
import io
import json
import os
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location(
    "fetch_inflammatory_arrays", ROOT / "scripts/reanalysis/fetch_inflammatory_arrays.py"
)
DOWNLOADER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DOWNLOADER)


class CompleteFileTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="pdac-download-test-")
        self.cache = Path(self.directory.name).resolve()
        self.assertEqual(self.cache.parent, Path(tempfile.gettempdir()).resolve())
        self.output = self.cache / "emtab1791_arrays"
        self.output.mkdir()
        with (ROOT / "data/publication/E-MTAB-1791.sdrf.txt").open(
            newline="", encoding="utf-8-sig"
        ) as stream:
            reader = csv.reader(stream, delimiter="\t")
            header, row = next(reader), next(reader)
        self.specimen = row[header.index("Source Name")]
        self.filename = row[header.index("Derived Array Data File")]
        self.target = self.output / self.filename
        self.marker = self.target.with_name(self.target.name + ".verified.json")
        self.payload = b"ProbeID\texpression\n1\t1.25\n"
        archive = io.BytesIO()
        with zipfile.ZipFile(archive, "w") as source:
            source.writestr(self.filename, self.payload)
        self.archive = archive.getvalue()

    def tearDown(self):
        self.directory.cleanup()

    def run_download(self, factory):
        argv = ["fetch", "--output-dir", str(self.output), self.specimen]
        with patch.dict(os.environ, {"PDAC_DATA_CACHE": str(self.cache)}), \
                patch.object(DOWNLOADER.sys, "argv", argv), \
                patch.object(DOWNLOADER, "RemoteZip", factory):
            DOWNLOADER.main()

    def archive_factory(self, url, **kwargs):
        return zipfile.ZipFile(io.BytesIO(self.archive))

    def test_success_publishes_verified_file(self):
        self.run_download(self.archive_factory)
        record = json.loads(self.marker.read_text(encoding="utf-8"))
        self.assertEqual(self.target.read_bytes(), self.payload)
        self.assertEqual(record["MD5"], hashlib.md5(self.payload).hexdigest())
        self.assertEqual(record["Bytes"], len(self.payload))
        self.assertFalse(self.target.with_name(self.target.name + ".part").exists())

    def test_full_size_interruption_is_not_published_and_retry_recovers(self):
        self.target.write_bytes(b"previous")
        self.marker.write_text("{}", encoding="utf-8")
        payload = self.payload

        class InterruptedStream(io.BytesIO):
            def read(self, size=-1):
                if self.tell() == len(payload):
                    raise OSError("Simulated transfer interruption before completion")
                return super().read(size)

        class InterruptedArchive(zipfile.ZipFile):
            def open(self, *args, **kwargs):
                return InterruptedStream(payload)

        with self.assertRaises(OSError):
            self.run_download(lambda url, **kwargs: InterruptedArchive(io.BytesIO(self.archive)))
        self.assertEqual(self.target.read_bytes(), b"previous")
        self.assertFalse(self.marker.exists())
        self.assertEqual(self.target.with_name(self.target.name + ".part").stat().st_size, len(payload))
        self.run_download(self.archive_factory)
        self.assertEqual(self.target.read_bytes(), payload)
        self.assertTrue(self.marker.exists())

    def test_bad_crc_is_not_published(self):
        corrupt = self.archive.replace(self.payload, b"X" + self.payload[1:], 1)
        self.assertNotEqual(corrupt, self.archive)
        with self.assertRaises(zipfile.BadZipFile):
            self.run_download(lambda url, **kwargs: zipfile.ZipFile(io.BytesIO(corrupt)))
        self.assertFalse(self.target.exists())
        self.assertFalse(self.marker.exists())


if __name__ == "__main__":
    unittest.main()
