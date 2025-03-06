# Automated test workflows

This directory stores workflow scripts that are picked up by GitHub Actions to automatically perform test installations of V-pipe on Mac OS and Linux systems and run end-to-end tests by executing tutorials with real example data.
This ensures successful installation and reproducible execution on different systems.
For each update of V-pipe these workflow scripts are automatically executed and report about installation problems or troubles on the test data

Some of these workflows rely on additional data present in the [`tests/` directory](../../tests).

## Tutorials

If you want to test your new V-pipe installation, we strongly encourage you to instead check our tutorials which provide real example data.

Tutorials for your first steps with V-pipe for different scenarios are available in the [docs/](../../docs/README.md) sub-directory.
