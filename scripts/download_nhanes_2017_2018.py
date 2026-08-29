from __future__ import annotations

from pathlib import Path
from urllib.request import Request, urlopen


BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles"
FILES = {
    "demo": "DEMO_J.xpt",
    "medical_conditions": "MCQ_J.xpt",
    "prescriptions": "RXQ_RX_J.xpt",
    "biochemistry": "BIOPRO_J.xpt",
    "body_measures": "BMX_J.xpt",
    "diabetes": "DIQ_J.xpt",
    "blood_pressure": "BPQ_J.xpt",
}


def _is_xport(payload: bytes) -> bool:
    # SAS transport v5/v8 files begin with an ASCII HEADER RECORD marker.
    return payload[:80].lstrip().startswith(b"HEADER RECORD")


def download(output_dir: str | Path = "data/nhanes_2017_2018") -> dict[str, Path]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for label, filename in FILES.items():
        url = f"{BASE}/{filename}"
        target = out / filename.upper()
        print(f"download {label}: {url}")
        request = Request(url, headers={"User-Agent": "rwe-programming-public-data-validator/1.0"})
        with urlopen(request, timeout=60) as response:
            payload = response.read()
            content_type = response.headers.get("Content-Type", "")
        if not _is_xport(payload):
            preview = payload[:120].decode("utf-8", errors="replace").replace("\n", " ")
            raise RuntimeError(
                f"CDC response for {filename} is not a SAS XPORT file "
                f"(content-type={content_type!r}, bytes={len(payload)}, preview={preview!r})"
            )
        target.write_bytes(payload)
        print(f"verified {filename}: {len(payload):,} bytes")
        paths[label] = target
    return paths


if __name__ == "__main__":
    download()
    print("NHANES public XPT files downloaded and signature-verified. Do not commit raw files.")
