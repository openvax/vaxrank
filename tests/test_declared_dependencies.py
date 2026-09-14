"""The installed environment must satisfy requirements.txt, and the floors
that exist for a stated reason must keep satisfying that reason.

Nothing else checks this. A dependency can drift below a declared floor and
the suite still runs, failing somewhere unrelated to the cause: mhcgnomes
3.33.4 surfaced as an ASCII-locale CLI failure inside a subprocess, and
mhctools 3.35.1 surfaced as a DSL parse error about a bracket. Both were
three imports away from the declared requirement they violated.
"""

from importlib.metadata import PackageNotFoundError, version
import importlib
from pathlib import Path

from packaging.requirements import Requirement
from packaging.version import Version
import pytest


REQUIREMENTS = Path(__file__).resolve().parent.parent / "requirements.txt"


def declared_requirements():
    """Parse requirements.txt into Requirement objects.

    Skips blank and comment lines, and anything not expressible as a plain
    name/specifier requirement (a VCS URL, an editable install), rather than
    guessing at what such a line demands.
    """
    parsed = []
    for raw in REQUIREMENTS.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        try:
            requirement = Requirement(line)
        except Exception:
            continue
        if requirement.url or requirement.marker:
            continue
        parsed.append(requirement)
    return parsed


def requirement_named(name):
    for requirement in declared_requirements():
        if requirement.name == name:
            return requirement
    raise AssertionError("requirements.txt declares no %s requirement" % name)


def floor_of(requirement):
    """Lowest version the specifier admits, or None if it sets no floor."""
    floors = [Version(spec.version) for spec in requirement.specifier
              if spec.operator in (">=", "==", "~=")]
    return max(floors) if floors else None


def test_requirements_txt_parses_into_at_least_the_known_dependencies():
    names = {requirement.name for requirement in declared_requirements()}
    # Guards the parser itself: a silently-empty parse would make every
    # other check in this module vacuously pass.
    assert {"isovar", "mhctools", "mhcgnomes", "topiary", "varcode"} <= names


#: Distribution name -> import name, where they differ.
_MODULE_NAMES = {"pyyaml": "yaml", "pillow": "PIL"}


def effective_version(distribution):
    """The version that actually governs behavior at runtime.

    An editable install's dist-info metadata records the version present
    when it was registered, which drifts from the source tree it points at.
    The imported module is what the suite runs against, so prefer its
    ``__version__`` and fall back to packaging metadata.
    """
    module_name = _MODULE_NAMES.get(distribution, distribution.replace("-", "_"))
    try:
        module = importlib.import_module(module_name)
    except Exception:
        module = None
    runtime = getattr(module, "__version__", None)
    if runtime:
        return str(runtime)
    try:
        return version(distribution)
    except PackageNotFoundError:
        return None


@pytest.mark.parametrize(
    "requirement", declared_requirements(), ids=lambda r: r.name)
def test_installed_distribution_satisfies_declared_requirement(requirement):
    installed = effective_version(requirement.name)
    if installed is None:
        pytest.skip("%s is not installed in this environment" % requirement.name)
    assert requirement.specifier.contains(installed, prereleases=True), (
        "installed %s %s violates the declared requirement %s; the suite would "
        "otherwise fail somewhere unrelated to this cause"
        % (requirement.name, installed, requirement))


def test_mhctools_floor_keeps_both_reasons_it_was_raised_for():
    """One pin, two independent reasons, neither enforced by its own comment.

    3.44.0 makes CleavageModel.scored_endpoint required for quantitative
    models, which vaxrank constructs; 3.39.0 is where mhctools' Kind class
    gained serum_half_life/blood_half_life, without which topiary's ranking
    DSL cannot reach either kind (openvax/topiary#300). Relaxing the floor
    for one reason must not silently drop below the other.
    """
    floor = floor_of(requirement_named("mhctools"))
    assert floor is not None, "mhctools must declare a lower bound"
    assert floor >= Version("3.44.0"), (
        "vaxrank constructs CleavageModel with scored_endpoint/motif_strictness, "
        "which only exist from mhctools 3.44.0")
    assert floor >= Version("3.39.0"), (
        "below mhctools 3.39.0 the ranking DSL cannot reach serum_half_life or "
        "blood_half_life (openvax/topiary#300)")


def test_mhctools_floor_matches_the_cleavage_contract_it_claims():
    """The floor's stated reason must be true of the installed mhctools."""
    from mhctools.cleavage import CleavageModel

    model = CleavageModel(
        "floor-check", "1", "proteasome", "", "", ("proteasomal",),
        "quantitative_model", ("declared-dependency-test",), "synthetic",
        "Not biological evidence", "probability", "dimensionless",
        scored_endpoint="site_cleavage")
    assert model.scored_endpoint == "site_cleavage"
    with pytest.raises(ValueError):
        CleavageModel(
            "floor-check", "1", "proteasome", "", "", ("proteasomal",),
            "quantitative_model", ("declared-dependency-test",), "synthetic",
            "Not biological evidence", "probability", "dimensionless")


def test_topiary_floor_keeps_the_cached_scan_fix_it_was_raised_for():
    """5.55.1 rebinds a cached protein scan to the requested occurrence.

    Below it, CachedPredictor can emit duplicate peptide_offset columns and
    crash the real-cache selection path (openvax/topiary#296).
    """
    floor = floor_of(requirement_named("topiary"))
    assert floor is not None and floor >= Version("5.55.1")
