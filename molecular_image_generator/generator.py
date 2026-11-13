"""Read molecular inputs and draw 2D PNG images."""

import csv
import colorsys
import hashlib
import json
import re
import tempfile
from pathlib import Path

from molecular_selector import FunctionalGroup, MolecularSelector
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from sdf_downloader.downloader import download_csv

INPUT_TYPES = ("sdf", "smiles", "name", "inchi")
_NON_GROUP_PATTERNS = frozenset({FunctionalGroup.ALKANE_CARBON, FunctionalGroup.ALPHA_AMINO_ACID})
_COMMON_AROMATIC_SMARTS = ("[a]", "[a]:[a]")
_GROUP_FAMILIES = {
    "amine": "amine",
    "primary amine": "amine",
    "secondary amine": "amine",
    "tertiary amine": "amine",
    "ammonium": "amine",
    "carboxylic acid": "carboxyl",
    "carboxylate": "carboxyl",
    "arene": "aromatic",
    "aromatic ring": "aromatic",
}
_DEFAULT_COLORS = {
    "amine": (0.35, 0.65, 1.0),
    "carboxyl": (1.0, 0.45, 0.35),
    "aromatic": (0.72, 0.55, 0.9),
}
_HIGHLIGHT_STRENGTH = 0.48


def _highlight_color(identity, attempt=0):
    if not attempt and identity in _DEFAULT_COLORS:
        return _DEFAULT_COLORS[identity]
    seed = identity if not attempt else f"{identity}:{attempt}"
    digest = hashlib.sha256(seed.encode("utf-8")).digest()
    hue = int.from_bytes(digest[:4], "big") / 2**32
    return colorsys.hsv_to_rgb(hue, 0.48 + digest[4] / 2550, 0.95)


def _csv_highlights(value, column):
    value = value.strip()
    if not value:
        return []
    after_bracket = value[1:].lstrip() if value.startswith("[") else ""
    is_array = value.startswith("[") and (
        column == "highlight" or after_bracket.startswith(('"', "]"))
    )
    if not is_array:
        return [value]
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{column} must contain a valid JSON array: {exc.msg}") from exc
    if not isinstance(parsed, list) or not all(isinstance(item, str) and item.strip() for item in parsed):
        raise ValueError(f"{column} must be a JSON array of nonempty strings")
    return parsed


def _compound_prefix(value):
    value = value.strip()
    if (value in (".", "..") or any(ord(char) < 32 or ord(char) == 127
                                    or char in '/\\:*?"<>|' for char in value)):
        raise ValueError("compound must be a filename without path separators or reserved characters")
    return value


class ImageGenerator:
    """Generate an individual PNG for each molecule in an input."""

    def __init__(self, *, width=600, height=600, highlight_names=(), highlight_smarts=(),
                 highlight_common=False):
        if type(width) is not int or type(height) is not int or width <= 0 or height <= 0:
            raise ValueError("width and height must be positive integers")
        if type(highlight_common) is not bool:
            raise ValueError("highlight_common must be a boolean")
        self.width = width
        self.height = height
        self.highlight_common = highlight_common
        self.selector = MolecularSelector()
        self.highlight_queries = []
        for key, values in (("name", highlight_names), ("smarts", highlight_smarts)):
            if isinstance(values, str):
                raise ValueError(f"highlight {key} values must be a sequence of strings")
            for value in values:
                if not isinstance(value, str) or not value.strip():
                    raise ValueError(f"highlight {key} must be a nonempty string")
                self.selector.select_molecules([Chem.MolFromSmiles("C")], **{key: value})
                self.highlight_queries.append((key, value))

    def _highlight_colors(self, mol):
        atom_groups, bond_groups = {}, {}
        queries = list(self.highlight_queries)
        if self.highlight_common:
            queries.extend(("name", group.value) for group in FunctionalGroup
                           if group not in _NON_GROUP_PATTERNS)
            queries.extend(("smarts", pattern) for pattern in _COMMON_AROMATIC_SMARTS)
        for key, value in queries:
            result = self.selector.select_molecules([mol], **{key: value})
            if key == "name":
                name = result["selector"]["name"]
                identity = _GROUP_FAMILIES.get(name, f"name:{name}")
            else:
                identity = "aromatic" if value in _COMMON_AROMATIC_SMARTS else f"smarts:{value}"
            for match in result["molecules"][0]["matches"]:
                matched = set(match["atom_indices"])
                for index in matched:
                    atom_groups.setdefault(index, set()).add(identity)
                for bond in mol.GetBonds():
                    if bond.GetBeginAtomIdx() in matched and bond.GetEndAtomIdx() in matched:
                        bond_groups.setdefault(bond.GetIdx(), set()).add(identity)
        identities = set().union(*atom_groups.values())
        colors, used = {}, set()
        for identity in sorted(identities, key=lambda value: (value not in _DEFAULT_COLORS, value)):
            attempt = 0
            color = _highlight_color(identity, attempt)
            rgb = tuple(int(channel * 255) for channel in color)
            while rgb in used:
                attempt += 1
                color = _highlight_color(identity, attempt)
                rgb = tuple(int(channel * 255) for channel in color)
            colors[identity] = color
            used.add(rgb)
        return tuple(
            {index: sorted(colors[identity] for identity in groups)
             for index, groups in group_map.items()}
            for group_map in (atom_groups, bond_groups)
        )

    def _highlight_indices(self, mol):
        atoms, bonds = self._highlight_colors(mol)
        return sorted(atoms), sorted(bonds)

    @staticmethod
    def _draw_highlighted(drawer, mol, atom_colors, bond_colors):
        options = drawer.drawOptions()
        options.highlightRadius = 0.22
        options.highlightBondWidthMultiplier = 5
        options.standardColoursForHighlightedAtoms = True

        def soften(color_map):
            return {
                index: [tuple(1 - (1 - channel) * _HIGHLIGHT_STRENGTH for channel in color)
                        for color in colors]
                for index, colors in color_map.items()
            }

        drawer.DrawMoleculeWithHighlights(mol, "", soften(atom_colors), soften(bond_colors), {}, {})

    @staticmethod
    def _molecules(kind, value):
        if kind == "sdf":
            yield from MolecularSelector.read_sdf(value)
        elif kind == "name":
            with tempfile.TemporaryDirectory(prefix="molecular-image-") as directory:
                root = Path(directory)
                csv_path = root / "input.csv"
                with csv_path.open("w", encoding="utf-8", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(["name"])
                    writer.writerow([value])
                output = root / "sdf"
                if download_csv(csv_path, output, record_type="2d"):
                    raise ValueError(f"could not resolve molecule name {value!r}")
                paths = list(output.glob("*.sdf"))
                if len(paths) != 1:
                    raise ValueError(f"expected one structure for molecule name {value!r}")
                yield from MolecularSelector.read_sdf(paths[0])
        elif kind in ("smiles", "inchi"):
            molecule = Chem.MolFromSmiles(value) if kind == "smiles" else Chem.MolFromInchi(value)
            if molecule is None or molecule.GetNumAtoms() == 0:
                raise ValueError(f"invalid {kind}: {value!r}")
            yield molecule
        else:
            raise ValueError(f"unsupported input type: {kind!r}")

    def generate(self, kind, value, output_dir, *, prefix="molecule"):
        """Return paths of PNG files created for a single input.

        SDF inputs may contain multiple records. Existing images are not overwritten.
        """
        if kind not in INPUT_TYPES:
            raise ValueError(f"unsupported input type: {kind!r}")
        if not isinstance(value, (str, Path)) or not str(value).strip():
            raise ValueError(f"{kind} input must not be empty")
        if not isinstance(prefix, str) or not prefix or Path(prefix).name != prefix or prefix in (".", ".."):
            raise ValueError("prefix must be a nonempty filename component")
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        paths = []
        try:
            for index, mol in enumerate(self._molecules(kind, str(value)), 1):
                if mol is None or mol.GetNumAtoms() == 0:
                    raise ValueError(f"{kind} record {index} contains no atoms")
                atoms, bonds = self._highlight_colors(mol)
                depiction = Chem.Mol(mol)
                rdDepictor.Compute2DCoords(depiction)
                drawer = rdMolDraw2D.MolDraw2DCairo(self.width, self.height)
                if atoms:
                    self._draw_highlighted(drawer, depiction, atoms, bonds)
                else:
                    drawer.DrawMolecule(depiction)
                drawer.FinishDrawing()
                path = output / f"{prefix}-{index}.png"
                with path.open("xb") as stream:
                    stream.write(drawer.GetDrawingText())
                paths.append(path)
        except (ValueError, OSError, RuntimeError) as exc:
            if paths:
                raise ValueError(f"{exc}; images already created: {', '.join(map(str, paths))}") from exc
            raise
        if not paths:
            raise ValueError(f"{kind} input contains no molecule records")
        return paths

    def generate_csv(self, csv_file, output_dir, *, delimiter=","):
        """Return (created paths, row errors), processing valid rows independently."""
        if len(delimiter) != 1:
            raise ValueError("CSV delimiter must be a single character")
        path = Path(csv_file)
        paths, errors = [], []
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(stream, delimiter=delimiter, strict=True)
            header = next(reader, None)
            if not header:
                raise ValueError("CSV must have a header")
            header = [column.strip().lower() for column in header]
            if len(set(header)) != len(header):
                raise ValueError("CSV contains duplicate headers")
            if "highlight_all" in header:
                raise ValueError("CSV column highlight_all was replaced by highlight_common")
            invalid_smarts_columns = [column for column in header
                                      if column.startswith("smarts_")
                                      and not re.fullmatch(r"smarts_[0-9]+", column)]
            if invalid_smarts_columns:
                raise ValueError(f"invalid numbered SMARTS columns: {', '.join(invalid_smarts_columns)}")
            smarts_columns = [column for column in header if re.fullmatch(r"smarts_[0-9]+", column)]
            if not set(header).intersection(INPUT_TYPES):
                raise ValueError(f"CSV must include at least one of: {', '.join(INPUT_TYPES)}")
            for row in reader:
                line = reader.line_num
                if not row:
                    continue
                try:
                    if len(row) != len(header):
                        raise ValueError(f"expected {len(header)} columns, got {len(row)}")
                    values = [(key, cell.strip()) for key, cell in zip(header, row)
                              if key in INPUT_TYPES and cell.strip()]
                    if len(values) != 1:
                        raise ValueError("row must contain exactly one input value")
                    kind, value = values[0]
                    columns = dict(zip(header, row))
                    compound = _compound_prefix(columns.get("compound", ""))
                    names = _csv_highlights(columns.get("highlight", ""), "highlight")
                    smarts = _csv_highlights(columns.get("highlight_smarts", ""), "highlight_smarts")
                    smarts.extend(columns[column].strip() for column in smarts_columns
                                  if columns[column].strip())
                    common_groups = columns.get("highlight_common", "").strip().lower()
                    if common_groups not in ("", "true", "false", "1", "0"):
                        raise ValueError("highlight_common must be true, false, 1, 0, or empty")
                    generator = self
                    if names or smarts or common_groups in ("true", "1"):
                        generator = ImageGenerator(
                            width=self.width,
                            height=self.height,
                            highlight_common=self.highlight_common or common_groups in ("true", "1"),
                            highlight_names=[v for k, v in self.highlight_queries if k == "name"] + names,
                            highlight_smarts=[v for k, v in self.highlight_queries if k == "smarts"] + smarts,
                        )
                    if kind == "sdf":
                        sdf_path = Path(value)
                        value = sdf_path if sdf_path.is_absolute() else path.resolve().parent / sdf_path
                    paths.extend(generator.generate(
                        kind, value, output_dir, prefix=compound or f"row-{line}"
                    ))
                except (ValueError, OSError, RuntimeError) as exc:
                    errors.append(f"CSV line {line}: {exc}")
        if not paths and not errors:
            raise ValueError("CSV contains no input rows")
        return paths, errors
