from __future__ import annotations

from pathlib import Path
from urllib.request import urlretrieve


BASE = "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018"
FILES = {
    "demo": "DEMO_J.XPT",
    "medical_conditions": "MCQ_J.XPT",
    "prescriptions": "RXQ_RX_J.XPT",
    "biochemistry": "BIOPRO_J.XPT",
    "body_measures": "BMX_J.XPT",
    "diabetes": "DIQ_J.XPT",
    "blood_pressure": "BPQ_J.XPT",
}


def main(output_dir: str = "data/nhanes_2017_2018") -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for label, filename in FILES.items():
        destination = out / filename
        if destination.exists():
            print(f"skip {label}: {destination}")
            continue
        url = f"{BASE}/{filename}"
        print(f"download {label}: {url}")
        urlretrieve(url, destination)
    print("NHANES public files downloaded. Do not commit raw XPT files to the repository.")


if __name__ == "__main__":
    main()
