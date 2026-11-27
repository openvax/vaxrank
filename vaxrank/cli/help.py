# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Readable argparse help and table styling via ``datacache.cli_help``."""

from datacache.cli_help import (
    ANSI_BOLD,
    ANSI_BOLD_CYAN,
    ANSI_CYAN,
    ANSI_DIM,
    ANSI_RESET,
    ANSI_YELLOW,
    DIM_STATUS_VALUES,
    OPTION_PATTERN,
    WARNING_STATUS_VALUES,
    HelpArgumentParser,
    HelpFormatter,
    classify_status,
    color_enabled,
    enable_argcomplete,
    style_help,
    style_status,
    style_table_header,
)

__all__ = [
    "ANSI_BOLD",
    "ANSI_BOLD_CYAN",
    "ANSI_CYAN",
    "ANSI_DIM",
    "ANSI_RESET",
    "ANSI_YELLOW",
    "DIM_STATUS_VALUES",
    "HelpArgumentParser",
    "HelpFormatter",
    "OPTION_PATTERN",
    "WARNING_STATUS_VALUES",
    "classify_status",
    "color_enabled",
    "enable_argcomplete",
    "style_help",
    "style_status",
    "style_table_header",
]
