#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
if [ "$#" -gt 1 ]; then
    echo "Usage: $0 [output-directory]" >&2
    exit 2
fi

venv="$root/.venv"
if [ ! -e "$venv" ] && [ ! -L "$venv" ]; then
    "$root/examples/setup-venv.sh" "$venv"
fi
if [ ! -f "$venv/bin/activate" ] || [ ! -x "$venv/bin/python" ]; then
    echo "Invalid virtual environment at $venv; inspect it before retrying." >&2
    exit 1
fi
. "$venv/bin/activate"
python="$venv/bin/python"

if [ "$#" -eq 1 ]; then
    mkdir -p -- "$1"
    output_dir=$(CDPATH= cd -- "$1" && pwd)
else
    output_dir=$(mktemp -d "${TMPDIR:-/tmp}/molecular-image-example.XXXXXXXX")
fi

cd "$root"
"$python" -m molecular_image_generator --smiles CCO -o "$output_dir/smiles"
"$python" -m molecular_image_generator --inchi 'InChI=1S/H2O/h1H2' -o "$output_dir/inchi"
"$python" -m molecular_image_generator --sdf examples/ethanethiol.sdf -o "$output_dir/sdf"
"$python" -m molecular_image_generator --csv examples/molecules.csv -o "$output_dir/csv"
"$python" -m molecular_image_generator --smiles 'CC(=O)O' \
    --highlight "carboxylic acid" -o "$output_dir/highlight-name"
"$python" -m molecular_image_generator --smiles 'CCO' \
    --highlight-smarts '[O]' -o "$output_dir/highlight-smarts"
for example in amino-acids saccharides nucleosides nitrogenous-bases glycerol-fatty-acids; do
    "$python" -m molecular_image_generator --csv "examples/$example.csv" \
        -o "$output_dir/$example"
done
printf 'Example output: %s\n' "$output_dir"
