"""Commands and timing of the runtime benchmark."""

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

import benchmark_direct as direct

CONFIG = dict(python="/venv/bin/python", rubric="/bin/rubric", dssr="/bin/dssr",
              rnaview="/bin/rnaview", bpnet="/bin/bpnet")


def test_programs_run_with_their_own_commands():
    work = Path("/tmp/run")
    cif = work / "TEST.cif"
    assert direct.command_for("dssr", cif, work, CONFIG) == \
        ["/bin/dssr", "-i=/tmp/run/TEST.cif", "--pair-only"]
    assert direct.command_for("rubric", cif, work, CONFIG) == \
        ["/bin/rubric", str(cif), "--no-download", "--cache-dir", str(work)]
    assert direct.command_for("barnaba", cif, work, CONFIG) == \
        ["/venv/bin/barnaba", "ANNOTATE", "--pdb", str(cif), "-o", "annotation"]
    for tool in direct.common.TOOLS:
        assert "run_one.py" not in " ".join(direct.command_for(tool, cif, work, CONFIG))


def test_input_unchanged_and_checks_outside_the_timed_run(tmp_path):
    (tmp_path / "logs").mkdir()
    source = tmp_path / "TEST.cif"
    original = b"data_original\n# Keep every model and chain exactly as supplied\n"
    source.write_bytes(original)
    events = []

    def launch(command, **kwargs):
        cif = Path(command[1])
        assert cif.is_symlink()
        assert cif.resolve() == source.resolve()
        assert cif.read_bytes() == original
        events.append("launch")

        def communicate(**kwargs):
            events.append("exit")
            return b'{"n_pairs":0,"pairs":[]}', b""
        return SimpleNamespace(returncode=0, communicate=communicate)

    def clock():
        events.append("clock")
        return 10.0 if len(events) == 1 else 12.5

    def validate(*args):
        events.append("validate")
        return []

    with patch.object(direct.subprocess, "Popen", side_effect=launch), \
         patch.object(direct.time, "perf_counter", side_effect=clock), \
         patch.object(direct, "validate_output", side_effect=validate):
        result = direct.run_trial(
            "rubric", dict(pdb="TEST", nucleotides=1, source_path=str(source)), 0,
            SimpleNamespace(timeout=10), {"rubric": "rubric", "python": "/venv/bin/python"},
            tmp_path, tmp_path, "originalhash")
    assert events == ["clock", "launch", "exit", "clock", "validate"]
    assert result["wall_seconds"] == 2.5
    assert result["status"] == "ok"
    assert source.read_bytes() == original


def test_recorded_commands_and_config_have_no_machine_paths():
    config = dict(python="/home/me/venv/bin/python", rubric="/home/me/build/rubric",
                  dssr="/opt/x/x3dna-dssr", rnaview="/home/me/rv/bin/rnaview",
                  bpnet="/home/me/bp/bpnet", packages={"numpy": "2.2.6"})
    work = Path("/tmp/rnad-abc")
    cif = work / "1ABC.cif"
    for tool in direct.common.TOOLS:
        shown = direct.public_command(direct.command_for(tool, cif, work, config), work, config)
        assert not any(part.startswith("/") for part in shown), (tool, shown)
    assert direct.public_command(direct.command_for("dssr", cif, work, config), work, config) == \
        ["x3dna-dssr", "-i=WORK/1ABC.cif", "--pair-only"]
    public = direct.public_config(config)
    assert not any(isinstance(v, str) and v.startswith("/") for v in public.values())
    assert public["packages"] == {"numpy": "2.2.6"}


def test_truncated_dssr_output_rejected():
    with pytest.raises(ValueError):
        direct.validate_output("dssr", Path("x.cif"), Path("."), "", "total number of base pairs: 2")
    assert direct.validate_output("dssr", Path("x.cif"), Path("."), "",
                                  "total number of base pairs: 0") == []
