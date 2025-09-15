# ExScope

ExScope is a command-line tool for visualizing read counts of a specific gene given its transcript.

## Installation

To install ExScope, clone this repository and run the following command from the project's root directory:

```bash
pip install -e .
```

## Usage

```bash
exscope --bam your_data.bam --gtf your_annotations.gff3 --transcript YOUR_TRANSCRIPT_ID -o output.png
```

## Changelog

### v0.2.0

*   Added a `--version` flag to the command-line tool.
*   The tool no longer hangs when no output file is specified. Instead, it prints a message to the console.
*   Updated the version number to 0.2.0.

### v0.1.0

*   Initial release.
