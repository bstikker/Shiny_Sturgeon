# Changelog

All notable changes to Shiny-Sturgeon are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.0.0] - 2026-11-21

First tagged release.

### Added
- Shiny GUI with a **Setup** tab (YAML config, saved per case to `<output_folder>/config.yml`) and a **Live Run** tab.
- Live monitoring of a MinKNOW `bam_pass/` folder for a chosen barcode, optionally including `unclassified/` reads.
- Per-iteration pipeline: merge and sort new BAMs, `modkit adjust-mods` and `modkit extract full`, cumulative `sturgeon inputtobed`, and classification.
- Sturgeon v1 classification with a confidence-over-time plot and a top-calls table.
- Optional Sturgeon v2 classification running side by side with v1, sharing the `inputtobed` step when both use the same probes file.
- Periodic and on-demand CNV plots (QDNAseq/DNAcopy) built from all reads processed so far.
- Resuming an interrupted run from an existing output folder.
- Configuration validation before a run starts (absolute paths, required files, v2 parameters).

### Known limitations
- No v1↔v2 class-name mapping, so v2 classes use fallback colors instead of the shared clinical color scheme.

[Unreleased]: https://github.com/bstikker/Shiny_Sturgeon/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/bstikker/Shiny_Sturgeon/releases/tag/v1.0.0
