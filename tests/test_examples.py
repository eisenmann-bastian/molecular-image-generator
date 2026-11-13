import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ExampleScriptTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name) / "project with spaces"
        shutil.copytree(ROOT / "examples", self.root / "examples")
        self.env = dict(os.environ)
        self.env["CALL_LOG"] = str(self.root / "calls.txt")

    def run_script(self, name, *args):
        return subprocess.run(
            [str(self.root / "examples" / name), *map(str, args)],
            cwd=self.root.parent, env=self.env, capture_output=True, text=True,
        )

    def test_runner_uses_existing_environment_and_all_inputs(self):
        bin_dir = self.root / ".venv" / "bin"
        bin_dir.mkdir(parents=True)
        (bin_dir / "activate").write_text("export EXAMPLE_VENV_ACTIVE=1\n", encoding="utf-8")
        python = bin_dir / "python"
        python.write_text(
            '#!/bin/sh\n'
            '[ "${EXAMPLE_VENV_ACTIVE:-}" = 1 ] || exit 9\n'
            'printf "%s\\n" "$*" >> "$CALL_LOG"\n',
            encoding="utf-8",
        )
        python.chmod(0o755)
        output = self.root.parent / "output with spaces"
        result = self.run_script("run-example.sh", output)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f"Example output: {output}", result.stdout)
        calls = (self.root / "calls.txt").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(calls), 11)
        for expected, line in zip(("--smiles CCO", "--inchi InChI=", "--sdf examples/",
                                   "--csv examples/molecules.csv"), calls):
            self.assertIn(expected, line)
        self.assertIn("--highlight carboxylic acid", calls[4])
        self.assertIn("--highlight-smarts [O]", calls[5])
        for line, name in zip(calls[6:], ("amino-acids", "saccharides", "nucleosides", "nitrogenous-bases",
                                         "glycerol-fatty-acids")):
            self.assertIn(f"--csv examples/{name}.csv", line)
            self.assertIn(f"-o {output / name}", line)

    def test_existing_or_invalid_environment_is_not_modified(self):
        env_dir = self.root / ".venv"
        env_dir.mkdir()
        sentinel = env_dir / "sentinel"
        sentinel.write_text("unchanged", encoding="utf-8")
        setup = self.run_script("setup-venv.sh")
        self.assertNotEqual(setup.returncode, 0)
        self.assertIn("already exists", setup.stderr)
        run = self.run_script("run-example.sh")
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("Invalid virtual environment", run.stderr)
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "unchanged")

    def test_runner_bootstraps_missing_environment(self):
        bootstrap = self.root / "bootstrap-python"
        bootstrap.write_text(
            "#!/bin/sh\n"
            'if [ "$1" = "-c" ]; then exit 0; fi\n'
            '[ "$1" = "-m" ] && [ "$2" = "venv" ] || exit 3\n'
            'mkdir -p "$3/bin"\n'
            'printf "export EXAMPLE_VENV_ACTIVE=1\\n" > "$3/bin/activate"\n'
            "cat > \"$3/bin/python\" <<'EOF'\n"
            "#!/bin/sh\n"
            'if [ "$1" != "-m" ] || [ "$2" != "pip" ]; then\n'
            '    [ "${EXAMPLE_VENV_ACTIVE:-}" = 1 ] || exit 9\n'
            "fi\n"
            'printf "%s\\n" "$*" >> "$CALL_LOG"\n'
            "EOF\n"
            'chmod +x "$3/bin/python"\n',
            encoding="utf-8",
        )
        bootstrap.chmod(0o755)
        self.env["PYTHON"] = str(bootstrap)
        result = self.run_script("run-example.sh", self.root.parent / "output")
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = (self.root / "calls.txt").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(calls), 12)
        self.assertEqual(calls[0], f"-m pip install -e {self.root}")
        self.assertIn("--csv examples/molecules.csv", calls[4])
        self.assertIn("--highlight-smarts [O]", calls[6])
        for line, name in zip(calls[7:], ("amino-acids", "saccharides", "nucleosides", "nitrogenous-bases",
                                         "glycerol-fatty-acids")):
            self.assertIn(f"--csv examples/{name}.csv", line)

    def test_usage_and_empty_environment_path(self):
        for name in ("run-example.sh", "setup-venv.sh"):
            result = self.run_script(name, "one", "two")
            self.assertEqual(result.returncode, 2)
            self.assertIn("Usage:", result.stderr)
        result = self.run_script("setup-venv.sh", "")
        self.assertEqual(result.returncode, 2)
        self.assertIn("must not be empty", result.stderr)

    def test_test_runner_uses_requested_python(self):
        python = self.root / "fake-python"
        python.write_text('#!/bin/sh\nprintf "%s\\n" "$*" > "$CALL_LOG"\n', encoding="utf-8")
        python.chmod(0o755)
        self.env["PYTHON"] = str(python)
        result = self.run_script("run-tests.sh", "-q")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            (self.root / "calls.txt").read_text(encoding="utf-8").strip(),
            "-m unittest discover -s tests -v -q",
        )


if __name__ == "__main__":
    unittest.main()
