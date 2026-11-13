import csv
import io
import os
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from molecular_selector import FunctionalGroup, MolecularSelector
from PIL import Image

from molecular_image_generator.generator import ImageGenerator, _DEFAULT_COLORS, _highlight_color
from molecular_image_generator import cli


class ImageGeneratorTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.output = self.root / "images"
        self.generator = ImageGenerator(width=240, height=180)

    def sdf(self, filename, smiles):
        path = self.root / filename
        writer = Chem.SDWriter(str(path))
        for value in smiles:
            writer.write(Chem.MolFromSmiles(value))
        writer.close()
        return path

    def assert_png(self, path):
        data = path.read_bytes()
        self.assertEqual(data[:8], b"\x89PNG\r\n\x1a\n")
        self.assertEqual(data[12:16], b"IHDR")
        self.assertEqual(struct.unpack(">II", data[16:24]), (240, 180))

    def test_biological_example_csvs_render_expected_structures(self):
        examples = Path(__file__).resolve().parents[1] / "examples"
        expected = {
            "amino-acids.csv": {
                "Glycine": "C2H5NO2", "Alanine": "C3H7NO2", "Valine": "C5H11NO2",
                "Leucine": "C6H13NO2", "Isoleucine": "C6H13NO2", "Serine": "C3H7NO3",
                "Threonine": "C4H9NO3", "Cysteine": "C3H7NO2S", "Methionine": "C5H11NO2S",
                "Aspartic acid": "C4H7NO4", "Glutamic acid": "C5H9NO4",
                "Asparagine": "C4H8N2O3", "Glutamine": "C5H10N2O3",
                "Lysine": "C6H14N2O2", "Arginine": "C6H14N4O2", "Histidine": "C6H9N3O2",
                "Phenylalanine": "C9H11NO2", "Tyrosine": "C9H11NO3",
                "Tryptophan": "C11H12N2O2", "Proline": "C5H9NO2",
            },
            "saccharides.csv": {
                "Glucose": "C6H12O6", "Fructose": "C6H12O6", "Ribose": "C5H10O5",
                "Deoxyribose": "C5H10O4", "Sucrose": "C12H22O11",
            },
            "nucleosides.csv": {
                "Adenosine": "C10H13N5O4", "Guanosine": "C10H13N5O5",
                "Cytidine": "C9H13N3O5", "Uridine": "C9H12N2O6", "Thymidine": "C10H14N2O5",
            },
            "nitrogenous-bases.csv": {
                "Adenine": "C5H5N5", "Guanine": "C5H5N5O", "Cytosine": "C4H5N3O",
                "Thymine": "C5H6N2O2", "Uracil": "C4H4N2O2",
            },
            "glycerol-fatty-acids.csv": {
                "Glycerol": "C3H8O3", "Lauric acid": "C12H24O2",
                "Myristic acid": "C14H28O2", "Palmitic acid": "C16H32O2",
                "Stearic acid": "C18H36O2", "Oleic acid": "C18H34O2",
                "Linoleic acid": "C18H32O2", "Alpha-linolenic acid": "C18H30O2",
            },
        }
        for filename, formulas in expected.items():
            with self.subTest(filename=filename):
                path = examples / filename
                with path.open(encoding="utf-8", newline="") as stream:
                    rows = list(csv.DictReader(stream))
                self.assertEqual(len(rows), len(formulas))
                self.assertEqual({row["compound"] for row in rows}, set(formulas))
                for row in rows:
                    with self.subTest(compound=row["compound"]):
                        molecule = Chem.MolFromSmiles(row["smiles"])
                        self.assertIsNotNone(molecule)
                        self.assertEqual(rdMolDescriptors.CalcMolFormula(molecule), formulas[row["compound"]])
                        highlighter = ImageGenerator(
                            highlight_names=[row["highlight"]] if row.get("highlight") else [],
                            highlight_smarts=[row["highlight_smarts"]] if row.get("highlight_smarts") else [],
                            highlight_common=row.get("highlight_common") == "true",
                        )
                        self.assertTrue(highlighter._highlight_indices(molecule)[0])
                        if row["compound"] in ("Oleic acid", "Linoleic acid", "Alpha-linolenic acid"):
                            double_bonds = [bond for bond in molecule.GetBonds()
                                            if bond.GetBondType() == Chem.BondType.DOUBLE
                                            and bond.GetStereo() != Chem.BondStereo.STEREONONE]
                            count = {"Oleic acid": 1, "Linoleic acid": 2, "Alpha-linolenic acid": 3}
                            self.assertEqual(len(double_bonds), count[row["compound"]])
                            self.assertTrue(all(bond.GetStereo() == Chem.BondStereo.STEREOZ
                                                for bond in double_bonds))
                paths, errors = self.generator.generate_csv(path, self.output / filename)
                self.assertEqual(errors, [])
                self.assertEqual(len(paths), len(formulas))
                self.assertEqual([image.name for image in paths],
                                 [f'{row["compound"]}-1.png' for row in rows])
                for image in paths:
                    self.assert_png(image)

    def test_csv_compound_names_fallback_and_multi_record(self):
        self.sdf("structures.sdf", ["CCO", "CO"])
        csv_path = self.root / "names.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["CoMpOuNd", "smiles", "sdf"])
            writer.writerow(["  L-Alanine β  ", "CCO", ""])
            writer.writerow(["", "CO", ""])
            writer.writerow(["Two records", "", "structures.sdf"])
        paths, errors = self.generator.generate_csv(csv_path, self.output)
        self.assertEqual(errors, [])
        self.assertEqual([path.name for path in paths],
                         ["L-Alanine β-1.png", "row-3-1.png", "Two records-1.png", "Two records-2.png"])
        for path in paths:
            self.assert_png(path)

    def test_csv_compound_invalid_or_duplicate_is_row_error(self):
        csv_path = self.root / "names.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["compound", "smiles"])
            writer.writerows([
                ["First", "C"], ["First", "O"], ["../outside", "C"], [r"..\outside", "C"],
                [".", "C"], ["bad:name", "C"], ["bad\nname", "C"], ["Next", "C"],
            ])
        paths, errors = self.generator.generate_csv(csv_path, self.output)
        self.assertEqual([path.name for path in paths], ["First-1.png", "Next-1.png"])
        self.assertEqual(len(errors), 6)
        self.assertIn("File exists", errors[0])
        self.assertTrue(all("compound must be" in error for error in errors[1:]))
        self.assertFalse((self.root / "outside-1.png").exists())

    def test_smiles_and_inchi_generate_png(self):
        for kind, value in (("smiles", "CCO"), ("inchi", "InChI=1S/H2O/h1H2")):
            with self.subTest(kind=kind):
                path, = self.generator.generate(kind, value, self.output, prefix=kind)
                self.assert_png(path)

    def test_named_and_smarts_highlight_full_matches_and_bonds(self):
        molecule = Chem.MolFromSmiles("CC(=O)OCCO")
        generator = ImageGenerator(highlight_names=["carboxylic acid", "ester"],
                                   highlight_smarts=["[O:1]"])
        atoms, bonds = generator._highlight_indices(molecule)
        self.assertEqual(atoms, [1, 2, 3, 4, 6])
        self.assertEqual(bonds, [1, 2, 3])
        self.assertEqual(generator._highlight_indices(Chem.MolFromSmiles("CC")), ([], []))

    def test_highlight_common_covers_multiple_groups_without_alkane_or_backbone(self):
        alanine = Chem.MolFromSmiles("NC(C)C(=O)O")
        serine = Chem.MolFromSmiles("NC(CO)C(=O)O")
        cysteine = Chem.MolFromSmiles("NC(CS)C(=O)O")
        generator = ImageGenerator(highlight_common=True)
        self.assertEqual(generator._highlight_indices(alanine), ([0, 3, 4, 5], [3, 4]))
        self.assertEqual(generator._highlight_indices(serine)[0], [0, 2, 3, 4, 5, 6])
        self.assertEqual(generator._highlight_indices(cysteine)[0], [0, 3, 4, 5, 6])
        self.assertEqual(generator._highlight_indices(Chem.MolFromSmiles("CCC")), ([], []))
        with self.assertRaisesRegex(ValueError, "boolean"):
            ImageGenerator(highlight_common="true")

    def test_highlight_common_equals_union_of_supported_group_matches(self):
        selector = MolecularSelector()
        generator = ImageGenerator(highlight_common=True)
        excluded = {FunctionalGroup.ALKANE_CARBON, FunctionalGroup.ALPHA_AMINO_ACID}
        for smiles in ("NC(CCSC)C(=O)O", "NC(CCCCN)C(=O)O",
                       "NC(Cc1ccccc1)C(=O)O", "NC(Cc1c[nH]cn1)C(=O)O"):
            with self.subTest(smiles=smiles):
                molecule = Chem.MolFromSmiles(smiles)
                atoms, bonds = set(), set()
                for group in FunctionalGroup:
                    if group in excluded:
                        continue
                    result = selector.select_molecules([molecule], name=group.value)
                    for match in result["molecules"][0]["matches"]:
                        matched = set(match["atom_indices"])
                        atoms.update(matched)
                        bonds.update(bond.GetIdx() for bond in molecule.GetBonds()
                                     if bond.GetBeginAtomIdx() in matched
                                     and bond.GetEndAtomIdx() in matched)
                for pattern in ("[a]", "[a]:[a]"):
                    for match in molecule.GetSubstructMatches(Chem.MolFromSmarts(pattern)):
                        atoms.update(match)
                        bonds.update(bond.GetIdx() for bond in molecule.GetBonds()
                                     if bond.GetBeginAtomIdx() in match and bond.GetEndAtomIdx() in match)
                self.assertEqual(generator._highlight_indices(molecule), (sorted(atoms), sorted(bonds)))

    def test_common_highlights_heteroaromatic_and_fused_rings(self):
        generator = ImageGenerator(highlight_common=True)
        for smiles in ("NC(Cc1c[nH]cn1)C(=O)O", "c1ccc2ccccc2c1",
                       "c1ccccc1-c2ccccc2"):
            with self.subTest(smiles=smiles):
                molecule = Chem.MolFromSmiles(smiles)
                atoms, bonds = generator._highlight_indices(molecule)
                aromatic_atoms = {atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetIsAromatic()}
                aromatic_bonds = {bond.GetIdx() for bond in molecule.GetBonds() if bond.GetIsAromatic()}
                self.assertTrue(aromatic_atoms.issubset(atoms))
                self.assertTrue(aromatic_bonds.issubset(bonds))
                if "-" in smiles:
                    self.assertTrue(all(bond.GetIdx() not in bonds for bond in molecule.GetBonds()
                                        if not bond.GetIsAromatic()))

    def test_all_amino_acid_examples_highlight_nitrogen_and_carboxyl(self):
        path = Path(__file__).resolve().parents[1] / "examples" / "amino-acids.csv"
        highlighter = ImageGenerator(highlight_common=True)
        backbone = Chem.MolFromSmarts("[N;X3][C;X4][C;X3](=[O])[O;H1]")
        with path.open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                with self.subTest(compound=row["compound"]):
                    self.assertEqual(row["highlight_common"], "true")
                    mol = Chem.MolFromSmiles(row["smiles"])
                    highlighted, _ = highlighter._highlight_indices(mol)
                    match = mol.GetSubstructMatch(backbone)
                    self.assertTrue(match)
                    self.assertTrue({match[0], *match[2:]}.issubset(highlighted))

    def test_csv_highlight_common_is_per_row_and_combines_with_global(self):
        csv_path = self.root / "all.csv"
        csv_path.write_text(
            "smiles,highlight_common,highlight_smarts\n"
            "CCO,TRUE,\nCCO,false,\nCCO,1,[C]\nCCO,nope,\nCCO,0,\n",
            encoding="utf-8",
        )
        generator = ImageGenerator(width=240, height=180)
        paths, errors = generator.generate_csv(csv_path, self.output)
        self.assertEqual([p.name for p in paths],
                         ["row-2-1.png", "row-3-1.png", "row-4-1.png", "row-6-1.png"])
        self.assertEqual(len(errors), 1)
        self.assertIn("CSV line 5: highlight_common", errors[0])
        self.assertNotEqual(paths[0].read_bytes(), paths[1].read_bytes())
        self.assertEqual(paths[1].read_bytes(), paths[3].read_bytes())
        self.assertNotEqual(paths[2].read_bytes(), paths[0].read_bytes())
        always = ImageGenerator(width=240, height=180, highlight_common=True)
        always_paths, errors = always.generate_csv(csv_path, self.root / "global")
        self.assertEqual(len(errors), 1)
        self.assertEqual(always_paths[1].read_bytes(), paths[0].read_bytes())

    def test_cli_highlight_common(self):
        with patch("molecular_image_generator.cli.configure_logging"), patch.object(cli, "logger"):
            status = cli.main(["--smiles", "CCO", "--highlight-common", "--width", "240",
                               "--height", "180", "-o", str(self.output)])
        self.assertEqual(status, 0)
        expected = ImageGenerator(width=240, height=180, highlight_common=True)
        reference, = expected.generate("smiles", "CCO", self.output, prefix="reference")
        self.assertEqual((self.output / "molecule-1.png").read_bytes(), reference.read_bytes())
        with patch("molecular_image_generator.cli.configure_logging"), patch.object(cli, "logger") as logger:
            status = cli.main(["--list-selectors", "--highlight-common"])
        self.assertEqual(status, 1)
        logger.error.assert_called_once()

    def test_csv_common_adds_to_global_highlights_for_multiple_sdf_records(self):
        self.sdf("structures.sdf", ["CCO", "CCO"])
        csv_path = self.root / "structures.csv"
        csv_path.write_text(
            "sdf,highlight_common,highlight_smarts\nstructures.sdf,true,[C]\n",
            encoding="utf-8",
        )
        generator = ImageGenerator(width=240, height=180, highlight_smarts=["[O]"])
        paths, errors = generator.generate_csv(csv_path, self.output)
        self.assertEqual(errors, [])
        self.assertEqual(len(paths), 2)
        reference_generator = ImageGenerator(width=240, height=180, highlight_common=True,
                                             highlight_smarts=["[O]", "[C]"])
        reference, = reference_generator.generate("smiles", "CCO", self.root / "reference")
        for path in paths:
            self.assertEqual(path.read_bytes(), reference.read_bytes())

        global_common = ImageGenerator(width=240, height=180, highlight_common=True)
        global_paths, errors = global_common.generate_csv(csv_path, self.root / "global")
        self.assertEqual(errors, [])
        reference_without_oxygen = ImageGenerator(width=240, height=180, highlight_common=True,
                                                  highlight_smarts=["[C]"])
        global_reference, = reference_without_oxygen.generate("smiles", "CCO", self.root / "global-ref")
        for path in global_paths:
            self.assertEqual(path.read_bytes(), global_reference.read_bytes())

    def test_common_group_colors_are_consistent_and_visible_in_png(self):
        generator = ImageGenerator(width=400, height=300, highlight_common=True)
        for smiles in ("NC(CO)C(=O)O", "NCC(=O)O", "O=C(O)C1CCCN1"):
            mol = Chem.MolFromSmiles(smiles)
            backbone = Chem.MolFromSmarts("[N;X3][C;X4][C;X3](=[O])[O;H1]")
            nitrogen, _, carbon, oxygen1, oxygen2 = mol.GetSubstructMatch(backbone)
            colors, bonds = generator._highlight_colors(mol)
            self.assertIn(_DEFAULT_COLORS["amine"], colors[nitrogen])
            for index in (carbon, oxygen1, oxygen2):
                self.assertIn(_DEFAULT_COLORS["carboxyl"], colors[index])
            self.assertIn(_DEFAULT_COLORS["carboxyl"],
                          bonds[mol.GetBondBetweenAtoms(carbon, oxygen1).GetIdx()])
        image_path, = generator.generate("smiles", "NC(CO)C(=O)O", self.output)
        with Image.open(image_path) as image:
            pixels = image.convert("RGB").getcolors(image.width * image.height)
            pixel_colors = {color for _, color in pixels}
        for color in (_DEFAULT_COLORS["amine"], _DEFAULT_COLORS["carboxyl"]):
            self.assertIn(tuple(round(channel * 255) for channel in color), pixel_colors)
        serine_colors, _ = generator._highlight_colors(Chem.MolFromSmiles("NC(CO)C(=O)O"))
        self.assertNotIn(serine_colors[2][0], (_DEFAULT_COLORS["amine"], _DEFAULT_COLORS["carboxyl"]))

    def test_overlapping_groups_retain_distinct_colors_and_deduplicate(self):
        molecule = Chem.MolFromSmiles("CC(=O)O")
        generator = ImageGenerator(highlight_names=["carboxylic acid", "carbonyl", "carboxylic acid"])
        atoms, bonds = generator._highlight_colors(molecule)
        self.assertEqual(len(atoms[1]), 2)
        self.assertEqual(len(atoms[2]), 2)
        self.assertEqual(len(bonds[1]), 2)
        self.assertIn(_DEFAULT_COLORS["carboxyl"], atoms[1])
        reversed_generator = ImageGenerator(highlight_names=["carbonyl", "carboxylic acid"])
        self.assertEqual((atoms, bonds), reversed_generator._highlight_colors(molecule))
        svg_drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
        depiction = Chem.Mol(molecule)
        rdDepictor.Compute2DCoords(depiction)
        svg_drawer.DrawMoleculeWithHighlights(depiction, "", atoms, bonds, {}, {})
        svg_drawer.FinishDrawing()
        svg = svg_drawer.GetDrawingText().upper()
        for color in atoms[1]:
            hex_color = "#{:02X}{:02X}{:02X}".format(*(int(channel * 255) for channel in color))
            self.assertIn(hex_color, svg)

    def test_same_group_aliases_and_multiple_matches_have_same_color(self):
        molecule = Chem.MolFromSmiles("NCCCCN")
        for name in ("amine", "primary amine", "primary_amine"):
            atoms, _ = ImageGenerator(highlight_names=[name])._highlight_colors(molecule)
            self.assertEqual(atoms[0], [_DEFAULT_COLORS["amine"]])
            self.assertEqual(atoms[5], [_DEFAULT_COLORS["amine"]])
        aromatic = ImageGenerator(highlight_common=True, highlight_names=["aromatic ring"])
        atoms, _ = aromatic._highlight_colors(Chem.MolFromSmiles("c1ccccc1"))
        self.assertTrue(all(colors == [_DEFAULT_COLORS["aromatic"]] for colors in atoms.values()))

    def test_color_collisions_are_resolved_without_changing_default_colors(self):
        def colliding_color(identity, attempt=0):
            if identity == "name:carbonyl" and not attempt:
                return _DEFAULT_COLORS["carboxyl"]
            return _highlight_color(identity, attempt)

        molecule = Chem.MolFromSmiles("CC(=O)O")
        with patch("molecular_image_generator.generator._highlight_color", side_effect=colliding_color):
            colors, _ = ImageGenerator(highlight_names=["carbonyl", "carboxylic acid"])._highlight_colors(molecule)
        self.assertEqual(len(colors[1]), 2)
        self.assertIn(_DEFAULT_COLORS["carboxyl"], colors[1])

    def test_numbered_smarts_are_additive_sparse_and_isolated_per_row(self):
        sdf = self.sdf("numbered.sdf", ["CCO", "CCO"])
        csv_path = self.root / "numbered.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["smiles", "sdf", "highlight_common", "highlight_smarts",
                             "SMARTS_0", "smarts_2", "smarts_42"])
            writer.writerow(["CCO", "", "false", '["[O]"]', "[C;H3]", "", "[C;H2]"])
            writer.writerow(["CCO", "", "", "", "", "", ""])
            writer.writerow(["CCO", "", "", "", "[", "", ""])
            writer.writerow(["", sdf.name, "true", "", "[C]", "", ""])
        generator = ImageGenerator(width=240, height=180)
        paths, errors = generator.generate_csv(csv_path, self.output)
        self.assertEqual([path.name for path in paths],
                         ["row-2-1.png", "row-3-1.png", "row-5-1.png", "row-5-2.png"])
        self.assertEqual(len(errors), 1)
        self.assertIn("CSV line 4", errors[0])
        self.assertIn("invalid SMARTS", errors[0])
        reference = ImageGenerator(width=240, height=180,
                                   highlight_smarts=["[O]", "[C;H3]", "[C;H2]"])
        expected, = reference.generate("smiles", "CCO", self.root / "expected")
        self.assertEqual(paths[0].read_bytes(), expected.read_bytes())
        plain, = self.generator.generate("smiles", "CCO", self.root / "plain")
        self.assertEqual(paths[1].read_bytes(), plain.read_bytes())
        common = ImageGenerator(width=240, height=180, highlight_common=True,
                                highlight_smarts=["[C]"])
        expected_common, = common.generate("smiles", "CCO", self.root / "common")
        self.assertEqual(paths[2].read_bytes(), expected_common.read_bytes())
        self.assertEqual(paths[3].read_bytes(), expected_common.read_bytes())

    def test_removed_option_and_malformed_smarts_columns_are_rejected(self):
        with self.assertRaises(TypeError):
            ImageGenerator(highlight_all=True)
        for header, error in (("highlight_all", "replaced by highlight_common"),
                              ("smarts_x", "invalid numbered SMARTS columns"),
                              ("smarts_-1", "invalid numbered SMARTS columns")):
            with self.subTest(header=header):
                csv_path = self.root / "invalid-header.csv"
                csv_path.write_text(f"smiles,{header}\nCCO,true\n", encoding="utf-8")
                with self.assertRaisesRegex(ValueError, error):
                    self.generator.generate_csv(csv_path, self.output)
        with patch("sys.stderr", new_callable=io.StringIO), self.assertRaises(SystemExit):
            cli.main(["--smiles", "CCO", "--highlight-all"])

    def test_aromatic_ring_and_multiple_separate_matches(self):
        generator = ImageGenerator(highlight_names=["aromatic ring"], highlight_smarts=["[OX2H]"])
        molecule = Chem.MolFromSmiles("Oc1ccccc1O")
        atoms, bonds = generator._highlight_indices(molecule)
        self.assertEqual(atoms, list(range(8)))
        self.assertEqual(len(bonds), 6)

    def test_highlight_changes_rendered_png(self):
        plain, = self.generator.generate("smiles", "CCO", self.output, prefix="plain")
        highlighted = ImageGenerator(width=240, height=180, highlight_names=["primary alcohol"])
        coloured, = highlighted.generate("smiles", "CCO", self.output, prefix="highlight")
        self.assert_png(coloured)
        self.assertNotEqual(plain.read_bytes(), coloured.read_bytes())
        no_match = ImageGenerator(width=240, height=180, highlight_names=["thiol"])
        unmatched, = no_match.generate("smiles", "CCO", self.output, prefix="unmatched")
        self.assertEqual(plain.read_bytes(), unmatched.read_bytes())

    def test_invalid_highlights_fail_before_output(self):
        for options in ({"highlight_names": ["not a selector"]}, {"highlight_smarts": ["["]},
                        {"highlight_smarts": [""]}, {"highlight_names": "thiol"}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                ImageGenerator(**options)
        self.assertFalse(self.output.exists())

    def test_highlights_apply_to_csv_and_multirecord_sdf(self):
        sdf = self.sdf("structures.sdf", ["CCO", "CCS"])
        csv_path = self.root / "input.csv"
        csv_path.write_text("sdf\nstructures.sdf\n", encoding="utf-8")
        generator = ImageGenerator(width=240, height=180, highlight_names=["thiol"],
                                   highlight_smarts=["[OX2H]"])
        paths, errors = generator.generate_csv(csv_path, self.output)
        self.assertEqual(errors, [])
        self.assertEqual([path.name for path in paths], ["row-2-1.png", "row-2-2.png"])
        reference = ImageGenerator(width=240, height=180)
        unhighlighted = reference.generate("sdf", sdf, self.output, prefix="plain")
        for coloured, plain in zip(paths, unhighlighted):
            self.assert_png(coloured)
            self.assertNotEqual(coloured.read_bytes(), plain.read_bytes())

    def test_csv_row_highlights_combine_with_global_and_do_not_leak(self):
        csv_path = self.root / "highlights.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["SmIlEs", "HIGHLIGHT", "highlight_smarts"])
            writer.writerow(["CCO", '["primary alcohol","ether"]', '["[O]","[C;H3]"]'])
            writer.writerow(["CCO", "", ""])
            writer.writerow(["CCO", "", "[#6,#8]"])
        generator = ImageGenerator(width=240, height=180, highlight_smarts=["[C;H2]"])
        paths, errors = generator.generate_csv(csv_path, self.output)
        self.assertEqual(errors, [])
        self.assertEqual([p.name for p in paths], ["row-2-1.png", "row-3-1.png", "row-4-1.png"])
        combined = ImageGenerator(
            width=240, height=180, highlight_names=["primary alcohol", "ether"],
            highlight_smarts=["[C;H2]", "[O]", "[C;H3]"],
        )
        reference, = combined.generate("smiles", "CCO", self.output, prefix="combined")
        global_only, = generator.generate("smiles", "CCO", self.output, prefix="global")
        self.assertEqual(paths[0].read_bytes(), reference.read_bytes())
        self.assertEqual(paths[1].read_bytes(), global_only.read_bytes())
        self.assertNotEqual(paths[2].read_bytes(), global_only.read_bytes())

    def test_csv_invalid_highlights_report_row_and_continue(self):
        csv_path = self.root / "invalid-highlights.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["smiles", "highlight", "highlight_smarts"])
            writer.writerow(["CCO", "unknown", ""])
            writer.writerow(["CCO", "", "["])
            writer.writerow(["CCO", '["thiol",', ""])
            writer.writerow(["CCO", '["thiol", 1]', ""])
            writer.writerow(["CCO", "", '["[O]", ""]'])
            writer.writerow(["CCO", "", "[O]"])
        paths, errors = self.generator.generate_csv(csv_path, self.output)
        self.assertEqual([path.name for path in paths], ["row-7-1.png"])
        self.assertEqual(len(errors), 5)
        for line, error in enumerate(errors, 2):
            self.assertIn(f"CSV line {line}:", error)
            self.assertFalse((self.output / f"row-{line}-1.png").exists())
        self.assertIn("JSON array", errors[2])
        self.assertIn("nonempty strings", errors[3])
        self.assert_png(paths[0])

    def test_csv_highlights_apply_to_each_sdf_record(self):
        self.sdf("structures.sdf", ["CCO", "CO"])
        csv_path = self.root / "input.csv"
        csv_path.write_text("sdf,highlight_smarts\nstructures.sdf,[O]\n", encoding="utf-8")
        paths, errors = self.generator.generate_csv(csv_path, self.output)
        self.assertEqual(errors, [])
        expected = ImageGenerator(width=240, height=180, highlight_smarts=["[O]"])
        refs = expected.generate("sdf", self.root / "structures.sdf", self.output, prefix="ref")
        self.assertEqual([p.read_bytes() for p in paths], [p.read_bytes() for p in refs])

    def test_cli_combines_csv_and_command_line_highlights(self):
        csv_path = self.root / "input.csv"
        csv_path.write_text("smiles,highlight_smarts\nCCO,[O]\n", encoding="utf-8")
        with patch("molecular_image_generator.cli.configure_logging"), patch.object(cli, "logger"):
            status = cli.main([
                "--csv", str(csv_path), "--highlight-smarts", "[C;H3]",
                "--width", "240", "--height", "180", "-o", str(self.output),
            ])
        self.assertEqual(status, 0)
        expected = ImageGenerator(width=240, height=180, highlight_smarts=["[C;H3]", "[O]"])
        reference, = expected.generate("smiles", "CCO", self.output, prefix="reference")
        self.assertEqual((self.output / "row-2-1.png").read_bytes(), reference.read_bytes())

    def test_sdf_multiple_records_and_invalid_input(self):
        path = self.sdf("molecules.sdf", ["CCO", "c1ccccc1"])
        images = self.generator.generate("sdf", path, self.output)
        self.assertEqual([p.name for p in images], ["molecule-1.png", "molecule-2.png"])
        for image in images:
            self.assert_png(image)
        with self.assertRaisesRegex(ValueError, "invalid smiles"):
            self.generator.generate("smiles", "not a molecule", self.output, prefix="invalid")
        with self.assertRaisesRegex(ValueError, "invalid inchi"):
            self.generator.generate("inchi", "not an inchi", self.output, prefix="invalid")
        with self.assertRaisesRegex(ValueError, "no molecule"):
            self.generator.generate("sdf", self.sdf("empty.sdf", []), self.output, prefix="empty")

    def test_name_uses_downloader_and_draws_locally(self):
        path = self.sdf("benzene.sdf", ["c1ccccc1"])

        def fake_download(csv_path, output, **kwargs):
            with Path(csv_path).open(newline="", encoding="utf-8") as source:
                self.assertEqual(list(csv.DictReader(source)), [{"name": "benzene"}])
            self.assertEqual(kwargs, {"record_type": "2d"})
            output.mkdir()
            (output / "row-2-cid-241.sdf").write_bytes(path.read_bytes())
            return 0

        with patch("molecular_image_generator.generator.download_csv", side_effect=fake_download):
            result = self.generator.generate("name", "benzene", self.output)
        self.assert_png(result[0])

    def test_name_failure_is_explicit(self):
        with patch("molecular_image_generator.generator.download_csv", return_value=1):
            with self.assertRaisesRegex(ValueError, "could not resolve molecule name"):
                self.generator.generate("name", "unknown", self.output)

    def test_name_with_real_downloader_and_mocked_pubchem(self):
        sdf_path = self.sdf("benzene.sdf", ["c1ccccc1"])
        text = sdf_path.read_text(encoding="utf-8").replace(
            "$$$$", "> <PUBCHEM_COMPOUND_CID>\n241\n\n$$$$"
        )
        with patch("sdf_downloader.downloader.pubchempy.get_sdf", return_value=text) as get_sdf:
            path, = self.generator.generate("name", "benzene", self.output)
        get_sdf.assert_called_once_with("benzene", namespace="name", record_type="2d")
        self.assert_png(path)

    def test_csv_rows_relative_sdf_and_errors(self):
        self.sdf("structures.sdf", ["C", "O"])
        csv_path = self.root / "input.csv"
        csv_path.write_text(
            "SMILES,sdf,name,inchi\nCCO,,,\n,structures.sdf,,\n"
            "bad,,,\nCC,O,,\n,,,InChI=1S/H2O/h1H2\n",
            encoding="utf-8",
        )
        images, errors = self.generator.generate_csv(csv_path, self.output)
        self.assertEqual([path.name for path in images],
                         ["row-2-1.png", "row-3-1.png", "row-3-2.png", "row-6-1.png"])
        self.assertEqual(len(errors), 2)
        self.assertIn("CSV line 4", errors[0])
        self.assertIn("exactly one", errors[1])
        for image in images:
            self.assert_png(image)

    def test_no_overwrite_and_validation(self):
        path, = self.generator.generate("smiles", "CC", self.output)
        original = path.read_bytes()
        with self.assertRaises(FileExistsError):
            self.generator.generate("smiles", "CO", self.output)
        self.assertEqual(path.read_bytes(), original)
        with self.assertRaises(ValueError):
            ImageGenerator(width=0)
        with self.assertRaises(ValueError):
            self.generator.generate("smiles", "", self.output)
        with self.assertRaises(ValueError):
            self.generator.generate("smiles", "CC", self.output, prefix="../bad")

    def test_multirecord_failure_reports_previously_created_images(self):
        sdf_path = self.sdf("structures.sdf", ["CC", "O"])
        sdf_path.write_bytes(sdf_path.read_bytes() + b"incomplete")
        with self.assertRaisesRegex(ValueError, "images already created:.*molecule-1.png.*molecule-2.png"):
            self.generator.generate("sdf", sdf_path, self.output)
        self.assert_png(self.output / "molecule-1.png")
        self.assert_png(self.output / "molecule-2.png")

    def test_cli_partial_failure_and_logging(self):
        path = self.root / "inputs.csv"
        path.write_text("smiles\nCC\nnot-a-smiles\n", encoding="utf-8")
        with patch("molecular_image_generator.cli.configure_logging"), patch.object(cli, "logger") as logger:
            status = cli.main(["--csv", str(path), "-o", str(self.output)])
        self.assertEqual(status, 1)
        self.assertTrue((self.output / "row-2-1.png").is_file())
        self.assertTrue(logger.error.called)

    def test_cli_accepts_multiple_highlights_and_lists_selectors(self):
        with patch("molecular_image_generator.cli.configure_logging"), patch.object(cli, "logger"):
            status = cli.main([
                "--smiles", "CCO", "--highlight", "primary alcohol",
                "--highlight-smarts", "[O]", "--highlight-smarts", "[C]",
                "--width", "240", "--height", "180", "-o", str(self.output),
            ])
        self.assertEqual(status, 0)
        self.assert_png(self.output / "molecule-1.png")
        with patch("molecular_image_generator.cli.configure_logging"), \
             patch("sys.stdout", new_callable=io.StringIO) as output:
            self.assertEqual(cli.main(["--list-selectors"]), 0)
        self.assertIn("primary alcohol", output.getvalue())

    def test_cli_rejects_unknown_selector_and_invalid_smarts(self):
        for option, value in (("--highlight", "unknown"), ("--highlight-smarts", "[")):
            with self.subTest(option=option), patch("molecular_image_generator.cli.configure_logging"), \
                 patch.object(cli, "logger") as logger:
                status = cli.main(["--smiles", "CCO", option, value, "-o", str(self.output)])
                self.assertEqual(status, 1)
                self.assertTrue(logger.error.called)
        self.assertFalse(self.output.exists())

    def test_module_cli_from_another_directory(self):
        project = Path(__file__).resolve().parents[1]
        env = dict(os.environ)
        inherited = [str(Path(part).resolve()) for part in env.get("PYTHONPATH", "").split(os.pathsep)
                     if part]
        env["PYTHONPATH"] = os.pathsep.join((str(project), *inherited))
        result = subprocess.run(
            [sys.executable, "-m", "molecular_image_generator", "--smiles", "CCO",
             "--width", "240", "--height", "180", "-o", str(self.output)],
            cwd=self.root, env=env, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assert_png(self.output / "molecule-1.png")
        self.assertIn("Saved", result.stderr)

    def test_csv_invalid_headers_and_empty(self):
        path = self.root / "input.csv"
        for content in ("x\n1\n", "sdf,SDF\nx,x\n", "smiles\n"):
            path.write_text(content, encoding="utf-8")
            with self.subTest(content=content), self.assertRaises(ValueError):
                self.generator.generate_csv(path, self.output)


if __name__ == "__main__":
    unittest.main()
