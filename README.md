# Molecular image generator

Generate one 2D PNG per molecular structure using RDKit. Images and their
2D coordinates are generated locally; PubChem image files are never used.

## Installation

Requires Python 3.10+, Git and the system dependencies required by RDKit.

```sh
python -m pip install .
```

The package links the GitHub versions of `molecular-selector`,
`molecular-utilities`, and `sdf-downloader`. **Before installing, publish the
packaging changes in the SDF downloader repository:** its existing GitHub
revision is not pip-installable. The other Git dependencies and RDKit are
installed by pip.

## Examples and test environment

Run the examples from any directory:

```sh
./examples/run-example.sh
./examples/run-example.sh /path/to/output-directory
```

The runner creates `.venv` and installs this checkout if the environment does
not exist, then activates it within the script. Existing environments are
used as-is; an incomplete one must be repaired or recreated manually. The
runner uses `.venv` even when another environment is active and does not
change the calling shell. Set `PYTHON` to choose the interpreter for creation.
The mixed CSV example looks up benzene on PubChem and requires network access.

The runner creates PNGs for SMILES, InChI, SDF and mixed CSV input, named
and SMARTS highlight examples, and all five biological CSV datasets below,
under separate output subdirectories. Without an argument it creates a
temporary output directory and prints its path. Existing images are not
overwritten.

The offline example CSVs contain the 20 standard amino acids,
representative saccharides, nucleosides, nitrogenous bases, and glycerol with
common fatty acids:

```sh
molecular-image-generator --csv examples/amino-acids.csv -o images-amino-acids
molecular-image-generator --csv examples/saccharides.csv -o images-saccharides
molecular-image-generator --csv examples/nucleosides.csv -o images-nucleosides
molecular-image-generator --csv examples/nitrogenous-bases.csv -o images-bases
molecular-image-generator --csv examples/glycerol-fatty-acids.csv -o images-lipids
```

The amino-acid, sugar and nucleoside examples do not specify stereochemistry;
they depict connectivity rather than a particular stereoisomer. The
unsaturated fatty-acid examples specify cis double bonds. The CSVs also
demonstrate per-row functional-group highlighting.

To create a new environment explicitly and run the offline tests:

```sh
./examples/setup-venv.sh                  # creates the project's .venv
# Or: ./examples/setup-venv.sh /path/to/new-environment
. .venv/bin/activate                     # or activate the chosen path
./examples/run-tests.sh
```

The setup script requires Python 3.10+, Git and network access to install the
GitHub dependencies. It refuses to modify an existing environment path. If
installation fails after creation, inspect or remove that environment before
retrying. The test runner uses the active interpreter or `PYTHON`.

## Usage

```sh
molecular-image-generator --smiles 'CCO' -o images
molecular-image-generator --inchi 'InChI=1S/H2O/h1H2' -o images-water
molecular-image-generator --name benzene -o images-benzene
molecular-image-generator --sdf molecules.sdf -o images-sdf
molecular-image-generator --csv molecules.csv -o images-batch
molecular-image-generator --csv examples/molecules.csv -o images-example
```

`--name` uses the SDF downloader to retrieve a single PubChem structure over
the network; lookup failures are reported. SMILES and InChI are parsed
locally. SDF files can contain multiple records; each receives an image.
`--width` and `--height` select pixel dimensions (default 600 each).

CSV headers are case-insensitive; supported columns are `sdf`, `smiles`,
`inchi` and `name`. Each row must contain **exactly one** nonempty supported
value. An optional `compound` column supplies the output filename prefix.
Optional `highlight`, `highlight_smarts`, `highlight_common`, and numbered
`smarts_0`, `smarts_1`, ... columns select per-row highlights. Other columns
are ignored. Relative SDF paths are relative to the CSV file.
For example:

```csv
name,smiles,inchi,sdf,highlight,highlight_smarts
benzene,,,,aromatic ring,
,CCO,,,primary alcohol,[O]
,,InChI=1S/H2O/h1H2,,,[O]
,,,molecules.sdf,thiol,
```

`--delimiter` changes the CSV separator. A standalone input creates
`molecule-1.png`, `molecule-2.png`, etc. CSV rows create
`<compound>-1.png`, `<compound>-2.png`, etc. when `compound` is nonempty,
or `row-<CSV line>-1.png`, `row-<CSV line>-2.png`, etc. otherwise. Spaces
and Unicode are retained; path separators, reserved filename characters and
control characters are rejected. Existing images are never overwritten, so
duplicate compound names or collisions with existing files cause row errors.
Invalid CSV rows are reported, other rows are processed,
and the command exits nonzero if any row failed. Successful images written
before an error remain in the output directory.
When an error occurs partway through a multi-record input, the error message
lists the images already created for that input.

The Python API provides `ImageGenerator.generate(kind, value, output_dir)`
and `generate_csv(csv_file, output_dir)`. The latter returns created paths and
row-error messages. Input molecule coordinates are not modified; a fresh 2D
depiction is generated on a copy.

Run offline tests with `python -m unittest discover -s tests -v`. Name lookups
are mocked in the tests.

## Highlight functional groups

Use named selectors from molecular-selector or provide SMARTS patterns:

```sh
molecular-image-generator --list-selectors
molecular-image-generator --smiles 'CC(=O)O' --highlight 'carboxylic acid' -o acid
molecular-image-generator --smiles 'CCO' --highlight-smarts '[O]' -o oxygen
molecular-image-generator --smiles 'NC(CO)C(=O)O' --highlight-common -o common-groups
molecular-image-generator --csv molecules.csv --highlight thiol --highlight-smarts '[O]' -o highlighted
```

The named and SMARTS options can be repeated and combined. Each detected
group type receives its own color, consistently across molecules; multiple
instances of the same type share a color. The default colors are blue for
amines and ammonium, coral for carboxylic acids and carboxylates, and purple
for aromatic atoms and bonds. Other named selectors and individual SMARTS
patterns receive deterministic colors. Identical colors are avoided within a
molecule, so in the rare case of a color collision a generated color can
depend on the other groups present. If groups overlap, all their colors are
shown on the shared atoms and bonds. All atoms of each matched group and
bonds between atoms within that same match are highlighted;
this is not limited to molecular-selector's narrower interaction-site atoms.
The options apply to every molecule, including every CSV row and SDF record.
No match produces an unhighlighted image. Unknown selectors and invalid
SMARTS are errors; command-line values are checked before processing inputs.

`--highlight-common` highlights every matching functional-group pattern in
molecular-selector except `alkane_carbon` (which matches ordinary saturated
carbon atoms) and `alpha_amino_acid` (which spans the amino-acid backbone).
It also highlights aromatic atoms and aromatic bonds using RDKit SMARTS
`[a]` and `[a]:[a]`, including heteroaromatic rings such as histidine's.
Matches and chemical coverage depend on RDKit's aromaticity perception and
molecular-selector's patterns;
arginine's side chain, for example, may be highlighted only in part. This
option does not imply chemically exhaustive functional-group detection.

CSV `highlight` and `highlight_smarts` cells may each contain one selector or
SMARTS string, or a JSON array of strings for multiple values. Quote JSON
arrays using standard CSV quoting. For example:

```csv
smiles,highlight,highlight_smarts
CCO,"[""primary alcohol"",""ether""]","[""[O]"",""[C;H3]""]"
```

For separate SMARTS columns, use as many `smarts_<number>` columns as needed.
Each nonempty cell is one SMARTS pattern; indices need not be consecutive.
These patterns are added to `highlight_smarts` and the other highlights.
CSV fields containing the delimiter must be quoted normally:

```csv
compound,smiles,highlight_common,smarts_0,smarts_1
Histidine,NC(Cc1c[nH]cn1)C(=O)O,true,[C;H2],[N]
```

Blank cells add nothing. Row highlights are added to command-line highlights;
they apply to every molecule in that row, including multi-record SDFs.
Invalid row highlights are reported for that row and subsequent rows still
run. A bracketed SMARTS such as `[O]` is a single SMARTS, not a JSON array.
The optional `highlight_common` column accepts `true`/`1` or `false`/`0`
(case-insensitive), or an empty cell. `true` enables common-group highlighting
for the row; `false` or empty does not disable a global `--highlight-common`.
The amino-acid example enables this for every row.
The old `highlight_all` CSV column, `--highlight-all` CLI option, and
`highlight_all` Python argument have been replaced by `highlight_common`.

For Python, pass sequences to the constructor, for example
`ImageGenerator(highlight_names=["thiol"], highlight_smarts=["[O]"],
highlight_common=True)`.
