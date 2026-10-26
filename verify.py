
"""Verify the retained release without training or modifying any files."""
from pathlib import Path
import hashlib
import json
ROOT=Path(__file__).resolve().parent

def verify():
    manifest=json.loads((ROOT/"published/manifest.json").read_text())
    for name,expected in manifest["files"].items():
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT): raise ValueError("Invalid manifest path")
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError("Retained artifact changed: "+name)
    import sys
    sys.path.insert(0,str(ROOT/"src"))
    import numpy as np
    from optics_portfolio.fiber import relative_error,sech
    with np.load(ROOT/"results/study_data.npz") as a:
        value=relative_error(a["field"],sech(a["t"])*np.exp(.5j*a["z"][:,None]))
    expected=json.loads((ROOT/"results/metrics.json").read_text())["soliton_relative_field_error"]
    np.testing.assert_allclose(value,expected,rtol=1e-10)
    assert json.loads((ROOT/"results/advanced_metrics.json").read_text())["certificate"]["passed"]
    return len(manifest["files"])

if __name__=="__main__": print("Verified",verify(),"retained artifact files; no retraining")
