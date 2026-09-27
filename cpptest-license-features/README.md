# C/C++test License Features

Queries the local or network license configured for one specific Parasoft
`cpptestcli` executable. It uses C/C++test's own fast `-check-license` operation
and does not analyze a project. The result explicitly identifies Automation,
Desktop CLI, both, or neither.

## Run

```powershell
python .\cpptest_license_features.py C:\parasoft\CPP_STD\cpptest\cpptestcli.exe
```

Use a separate settings file when the selected installation does not contain
the license configuration:

```powershell
python .\cpptest_license_features.py C:\path\cpptestcli.exe --settings C:\path\license.properties
```

For automation, add `--json`. Exit code `0` means the query succeeded. Exit
code `2` means the CLI could not be run, timed out, could not reach or acquire
its license, or returned an unrecognized response.

`Not recognized` means that the selected C/C++test version does not expose
that name through `-check-license`. It does not mean that a similarly named
License Server bundle or informational feature is unlicensed.

## Build a standalone Windows executable

```powershell
pyinstaller --onefile --name cpptest-license-features .\cpptest_license_features.py
.\dist\cpptest-license-features.exe C:\path\cpptestcli.exe
```

## Run the Linux executable

```bash
./dist/linux/cpptest-license-features \
  /opt/parasoft/cpptest/bin/cli/cpptestcli \
  --settings /path/to/license.properties
```

For C/C++test 2025.2, select `cpptest/bin/cli/cpptestcli`. The top-level
`cpptest/cpptestcli` launcher requires an analysis configuration and does not
provide the fast license-check operation.

The utility never reads or prints license passwords. The selected
`cpptestcli` loads its normal installation settings, plus `--settings` when
provided.
