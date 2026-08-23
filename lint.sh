#!/bin/bash
set -o errexit

ruff check vaxrank tests release_upload.py
echo 'Passes ruff check'
