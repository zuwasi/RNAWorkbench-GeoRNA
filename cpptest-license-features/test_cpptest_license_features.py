import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import cpptest_license_features as app


class ProbeTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.executable = Path(self.temp_dir.name) / "cpptestcli.exe"
        self.executable.touch()

    def tearDown(self):
        self.temp_dir.cleanup()

    @patch("cpptest_license_features.subprocess.run")
    def test_reports_partial_license(self, run):
        missing = set(app.FEATURES) - {"C++test", "MISRA C 2012 Rules", "SEI CERT C Rules"}
        run.return_value = subprocess.CompletedProcess(
            [],
            134,
            "Parasoft C/C++test Standard 2025.1\nLicense features not set: " + ", ".join(missing),
            "",
        )

        result = app.probe(self.executable)

        self.assertEqual(
            result.licensed,
            ("C++test", "MISRA C 2012 Rules", "SEI CERT C Rules"),
        )
        self.assertIn("MISRA C++ 2023 Rules", result.not_licensed)

    @patch("cpptest_license_features.subprocess.run")
    def test_removes_unknown_features_and_retries(self, run):
        run.side_effect = (
            subprocess.CompletedProcess(
                [], 134, "ERROR: Unknown license feature names in the -check-license option: MISRA C 2025 Rules", ""
            ),
            subprocess.CompletedProcess(
                [], 0, "Parasoft C/C++test Standard 2025.1\nAll license features are set: ok", ""
            ),
        )

        result = app.probe(self.executable)

        self.assertEqual(result.not_recognized, ("MISRA C 2025 Rules",))
        self.assertNotIn("MISRA C 2025 Rules", result.licensed)
        self.assertEqual(run.call_count, 2)

    @patch("cpptest_license_features.subprocess.run")
    def test_connection_failure_is_not_reported_as_no_features(self, run):
        run.return_value = subprocess.CompletedProcess(
            [],
            134,
            "License: Activation failed: Unable to connect to License Server\nLicense features not set: C++test",
            "",
        )

        with self.assertRaisesRegex(app.ProbeError, "Unable to connect"):
            app.probe(self.executable)

    @patch("cpptest_license_features.subprocess.run")
    def test_uses_separate_temporary_workspace(self, run):
        run.return_value = subprocess.CompletedProcess(
            [], 0, "Parasoft C/C++test Standard 2025.1\nAll license features are set: ok", ""
        )

        app.probe(self.executable)

        command = run.call_args.args[0]
        self.assertEqual(command[1], "-workspace")
        self.assertNotEqual(Path(command[2]), Path.cwd())


class LicenseModeTests(unittest.TestCase):
    def test_reports_both_when_both_are_licensed(self):
        self.assertEqual(
            app.classify_license_mode(("Automation", "Desktop Command Line")),
            "Both Automation and Desktop CLI",
        )

    def test_reports_automation_only(self):
        self.assertEqual(
            app.classify_license_mode(("Automation",)),
            "Automation edition",
        )

    def test_reports_desktop_cli_only(self):
        self.assertEqual(
            app.classify_license_mode(("Desktop Command Line",)),
            "Desktop CLI only",
        )

    def test_reports_neither(self):
        self.assertEqual(
            app.classify_license_mode(("Static Analysis",)),
            "Neither Automation nor Desktop CLI",
        )


if __name__ == "__main__":
    unittest.main()
