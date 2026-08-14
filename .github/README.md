# GitHub Configuration

This directory contains GitHub-specific configuration to ensure proper repository setup with security safeguards.

## .gitignore - Critical Security File

The `.gitignore` file is pre-configured to exclude:
- Input data files (*.fastq.gz, *.fasta, etc.)
- Configuration files with sensitive paths (cleangene.env)
- ARC infrastructure files
- Local development artifacts

## GitHub Workflows

GitHub Actions are configured for:
- Automated testing on push/pull request
- Python linting and formatting checks
- Dependency security scanning

All workflows use placeholder data and never access real sensitive information.

## Repository Setup Checklist

1. [ ] Create GitHub repository (public)
2. [ ] Set up branch protection rules
3. [ ] Configure CODEOWNERS file
4. [ ] Set up issue templates
5. [ ] Configure security advisories
6. [ ] Add CONTRIBUTING.md with data safety guidelines