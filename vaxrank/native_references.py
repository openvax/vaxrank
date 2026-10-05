"""Verified relative resources for native custom annotation references."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

from pyensembl import Genome

from .native_serialization import _json_loads


REFERENCE_FIELD = "__vaxrank_reference__"
REFERENCE_SCHEMA = "vaxrank.native_reference.v1"


def _digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class NativeReferences:
    """Bundle available files; restore identities without acquiring resources."""

    def __init__(self, path, previous=None):
        self.root = Path(path).resolve().parent
        self.directory = Path(Path(path).name + ".references")
        self.availability = []
        self.genomes = {}
        self._digests = {}
        self._expected = {item['resolved_path']: item for item in
                          (previous.availability if previous is not None else ())}

    @classmethod
    def combine(cls, references):
        """Retain reference identities when native sources are composed."""
        references = [reference for reference in references if reference is not None]
        if not references:
            return None
        combined = deepcopy(references[0])
        for reference in references[1:]:
            combined.availability.extend(reference.availability)
            combined.genomes.update(reference.genomes)
        return combined

    def _file_digest(self, path):
        stat = path.stat()
        key = (path, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
        if key not in self._digests:
            self._digests[key] = _digest(path)
        return self._digests[key]

    def _relative_path(self, relative):
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("Native reference path must stay inside the dataset directory")
        resolved = (self.root / path).resolve()
        if not resolved.is_relative_to(self.root):
            raise ValueError("Native reference path escapes the dataset directory")
        return resolved

    def _bundle_file(self, source):
        digest = self._file_digest(source)
        suffix = "".join(source.suffixes[-2:]) if source.suffix == ".gz" else source.suffix
        relative = self.directory / (digest + suffix)
        destination = self._relative_path(relative)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if self._file_digest(destination) != digest:
                raise ValueError(f"Native reference bundle is corrupt: {destination}")
        else:
            with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as stream:
                temporary = Path(stream.name)
            try:
                shutil.copyfile(source, temporary)
                if _digest(temporary) != digest:
                    raise ValueError(f"Reference changed during export: {source}")
                temporary.replace(destination)
            finally:
                temporary.unlink(missing_ok=True)
        return dict(path=str(relative), sha256=digest, size=destination.stat().st_size)

    def _resource(self, source):
        expected = self._expected.get(str(source.resolve()))
        if source.is_file():
            if expected is not None and (source.stat().st_size != expected['size']
                    or self._file_digest(source) != expected['sha256']):
                raise ValueError(f"Native reference checksum mismatch: {source}")
            return self._bundle_file(source)
        if expected is not None:
            # Evidence-only exports retain the identity of absent resources;
            # absence must not silently erase a previously recorded checksum.
            return dict(path=str(self.directory / Path(expected['path']).name),
                        sha256=expected['sha256'], size=expected['size'])
        return None

    def _encode_genome(self, value):
        metadata = value["__class__"]
        state = {key: item for key, item in value.items() if key != "__class__"}
        resources = []
        if metadata["__name__"] == "Genome":
            genome = Genome(**state)
            # Public inspection does not download, copy, index or open a DB.
            files = iter(genome.required_local_files())
            fields = []
            if state.get("genome_fasta_path_or_url"):
                fields.append(("genome_fasta_path_or_url", None))
            if state.get("gtf_path_or_url"):
                fields.append(("gtf_path_or_url", None))
            for key in ("transcript_fasta_paths_or_urls", "protein_fasta_paths_or_urls"):
                fields.extend((key, index) for index, _ in enumerate(state.get(key) or []))
            for field, index in fields:
                source = Path(next(files))
                resource = self._resource(source)
                if resource is not None:
                    resources.append(dict(field=field, index=index, **resource))
        else:
            # Exact Ensembl release/species is already portable. Attached local
            # reference DNA is the only source path in its public state.
            source = state.get("genome_fasta")
            resource = self._resource(Path(source)) if isinstance(source, str) else None
            if resource is not None:
                resources.append(dict(field="genome_fasta", index=None, **resource))
        if resources:
            value[REFERENCE_FIELD] = dict(schema=REFERENCE_SCHEMA, resources=resources)
        return value

    def _decode_genome(self, value):
        manifest = value.pop(REFERENCE_FIELD, None)
        if manifest is None:
            return value  # Legacy references retain their exact identity.
        if not isinstance(manifest, dict) or manifest.get("schema") != REFERENCE_SCHEMA:
            raise ValueError("Unsupported native reference manifest")
        resources = manifest.get("resources")
        if not isinstance(resources, list) or not resources:
            raise ValueError("Native reference manifest requires resources")
        custom = value["__class__"]["__name__"] == "Genome"
        seen = set()
        for resource in resources:
            if not isinstance(resource, dict) or set(resource) != {"field", "index", "path", "sha256", "size"}:
                raise ValueError("Malformed native reference resource")
            digest, size = resource["sha256"], resource["size"]
            if (not isinstance(digest, str) or len(digest) != 64
                    or any(c not in "0123456789abcdef" for c in digest)
                    or type(size) is not int or size < 0):
                raise ValueError("Invalid native reference checksum or size")
            field, index = resource["field"], resource["index"]
            if field not in ({"gtf_path_or_url", "genome_fasta_path_or_url",
                              "transcript_fasta_paths_or_urls", "protein_fasta_paths_or_urls"}
                             if custom else {"genome_fasta"}):
                raise ValueError("Invalid native reference resource field")
            key = (field, index)
            if key in seen:
                raise ValueError("Duplicate native reference resource")
            seen.add(key)
            path = self._relative_path(resource["path"])
            available = path.is_file()
            if available and (path.stat().st_size != resource["size"]
                              or self._file_digest(path) != resource["sha256"]):
                raise ValueError(f"Native reference checksum mismatch: {path}")
            self.availability.append(dict(resource, available=available, resolved_path=str(path)))
            if field.endswith("_paths_or_urls"):
                if type(index) is not int or not 0 <= index < len(value.get(field) or []):
                    raise ValueError("Invalid native reference resource index")
                value[field][index] = str(path)
            else:
                if index is not None:
                    raise ValueError("Invalid native reference scalar index")
                value[field] = str(path)
        if custom:
            # Derived indexes belong to this relocated bundle, never the old
            # annotation cache. Their creation remains an explicit operation.
            value["cache_directory_path"] = str(self.root / self.directory / "indexes")
            value["copy_local_files_to_cache"] = False
            value["decompress_on_download"] = False
            state = {key: item for key, item in value.items() if key != "__class__"}
            self.genomes[json.dumps(state, sort_keys=True)] = state
        return value

    def _walk(self, value, operation):
        if isinstance(value, list):
            return [self._walk(item, operation) for item in value]
        if not isinstance(value, dict):
            return value
        metadata = value.get("__class__", {})
        if isinstance(metadata, dict) and (metadata.get("__module__"), metadata.get("__name__")) in {
                ("pyensembl.genome", "Genome"), ("pyensembl.ensembl_release", "EnsemblRelease")}:
            value = operation(value)
        return {key: self._walk(item, operation) for key, item in value.items()}

    def transform(self, payload, *, decode=False):
        """Transform only native object fields, leaving source evidence intact."""
        operation = self._decode_genome if decode else self._encode_genome
        def native(value):
            return json.dumps(self._walk(_json_loads(value), operation)) if value is not None else None
        payload = deepcopy(payload)
        for field in ("epitopes", "provenance"):
            payload[field] = [native(value) for value in payload[field]]
        for field in ("antigens", "mutation_fragments"):
            payload[field] = {key: native(value) for key, value in payload.get(field, {}).items()}
        for source in payload.get("direct_sources", []):
            source["variant"] = native(source["variant"])
        for report in payload.get("construction_reports", []):
            report["genome"] = native(report.get("genome"))
            report["rows"] = native(report["rows"])
            report["records"] = [native(value) for value in report["records"]]
        return payload

    def index(self):
        """Prepare verified local annotations explicitly, without downloads."""
        missing = []
        for item in self.availability:
            path = Path(item["resolved_path"])
            if not path.is_file():
                missing.append(str(path))
            elif path.stat().st_size != item["size"] or self._file_digest(path) != item["sha256"]:
                raise ValueError(f"Native reference checksum mismatch: {path}")
        if missing:
            raise ValueError("Native reference resources are unavailable: " + ", ".join(sorted(set(missing))))
        for state in self.genomes.values():
            genome = Genome(**state)
            if not genome.required_local_files_exist():
                raise ValueError("Native annotation requires unavailable resources: " + str(genome))
            genome.index()
