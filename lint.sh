#!/bin/bash
set -o errexit

ruff check vaxrank tests
echo 'Passes ruff check'
